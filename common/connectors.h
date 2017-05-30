#pragma once

// Standard C++ includes
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <vector>

// GeNN includes
#include "sparseProjection.h"

//----------------------------------------------------------------------------
// SimpleSparseProjection
//----------------------------------------------------------------------------
class SimpleSparseProjection
{
private:
    typedef std::vector<unsigned int> Row;

public:
    SimpleSparseProjection()
    {
    }

    SimpleSparseProjection(unsigned int numRows) : m_RowIndices(numRows)
    {
    }

    //----------------------------------------------------------------------------
    // Public API
    //----------------------------------------------------------------------------
    unsigned int getNumRows() const
    {
        return m_RowIndices.size();
    }

    unsigned int calcNumSynapses() const
    {
        return std::accumulate(m_RowIndices.begin(), m_RowIndices.end(), 0,
                                [](unsigned int accumulator, const Row &row)
                                {
                                    return (accumulator + row.size());
                                });
    }

    unsigned int calcMaxColumns() const
    {
        // Find longest row
        auto longestRow = std::max_element(m_RowIndices.begin(), m_RowIndices.end(),
                                           [](const Row &a, const Row &b)
                                           {
                                               return (a.size() < b.size());
                                           });

        if(longestRow == m_RowIndices.end()) {
            return 0;
        }
        else {
            return longestRow->size();
        }

    }

    void setNumRows(unsigned int numRows)
    {
        m_RowIndices.resize(numRows);
    }

    void reserveColumns(unsigned int i, unsigned int numColumns)
    {
        checkRow(i);
        m_RowIndices[i].reserve(numColumns);
    }

    void addSynapse(unsigned int i, unsigned int j)
    {
        checkRow(i);
        m_RowIndices[i].push_back(j);
    }

    size_t rowLength(unsigned int i) const
    {
        checkRow(i);
        return m_RowIndices[i].size();
    }

    Row::const_iterator rowBegin(unsigned int i) const
    {
        checkRow(i);
        return m_RowIndices[i].cbegin();
    }

    Row::const_iterator rowEnd(unsigned int i) const
    {
        checkRow(i);
        return m_RowIndices[i].cend();
    }

    Row::iterator rowBegin(unsigned int i)
    {
        checkRow(i);
        return m_RowIndices[i].begin();
    }

    Row::iterator rowEnd(unsigned int i)
    {
        checkRow(i);
        return m_RowIndices[i].end();
    }

private:
    //----------------------------------------------------------------------------
    // Private API
    //----------------------------------------------------------------------------
    void checkRow(unsigned int i) const
    {
        if(i >= m_RowIndices.size()) {
            throw std::runtime_error("Invalid row");
        }
    }

    //----------------------------------------------------------------------------
    // Members
    //----------------------------------------------------------------------------
    std::vector<Row> m_RowIndices;
};

inline std::ostream & operator << (std::ostream &os, const SimpleSparseProjection& p)
{
    for(unsigned int i = 0; i < p.getNumRows(); i++)
    {
        os << i << ":";
        for(auto j = p.rowBegin(i); j != p.rowEnd(i); ++j)
        {
            os << *j << ",";
        }
        os << std::endl;
    }

    return os;
}

//----------------------------------------------------------------------------
// Typedefines
//----------------------------------------------------------------------------
typedef void (*AllocateFn)(unsigned int);

//----------------------------------------------------------------------------
// Functions
//----------------------------------------------------------------------------
void convertToSparseProjection(const SimpleSparseProjection &projection,
                               SparseProjection &sparseProjection, AllocateFn allocate)
{
    // Allocate projection
    allocate(projection.calcNumSynapses());

    // Loop through rows
    unsigned int s = 0;
    for(unsigned int i = 0; i < projection.getNumRows(); i++)
    {
        // Set synapse index of start of row
        sparseProjection.indInG[i] = s;

        // Copy row synapses into sparse projection
        std::copy(projection.rowBegin(i), projection.rowEnd(i),
                  &sparseProjection.ind[s]);

        // Add row length to total number of synapses
        s += projection.rowLength(i);
    }

    // Set final synapse index for end of last row
    sparseProjection.indInG[projection.getNumRows()] = s;
}

template<unsigned int N>
void convertToPartitionedSparseProjection(const SimpleSparseProjection &projection, unsigned int partitionSize,
                                          PartitionedSparseProjection<N> &sparseProjection, AllocateFn allocate)
{
    // Allocate projection
    allocate(projection.calcNumSynapses());

    // Loop through rows
    unsigned int s = 0;
    for(unsigned int i = 0; i < projection.getNumRows(); i++)
    {
        // Copy row synapses into sparse projection
        unsigned int *rowBegin = &sparseProjection.postIndices[s];
        unsigned int *rowEnd = std::copy(projection.rowBegin(i), projection.rowEnd(i), rowBegin);

        // Loop through partitions
        unsigned int *subRowPartitionBegin = rowBegin;
        for(unsigned int p = 0; p < N; p++)
        {
            // Calculate pivot point for this partition
            const unsigned int pivot = (p + 1) * partitionSize;

            // Partition remaining synapses in subrow so ones before pivot are moved to front
            unsigned int *subRowPartitionEnd = std::partition(subRowPartitionBegin, rowEnd,
                                                              [pivot](unsigned int j)
                                                              {
                                                                  return (j < pivot);

                                                              });

            // Add start and length indices to subrow
            sparseProjection.subRowBeginIndices[p][i] = s + (subRowPartitionBegin - rowBegin);
            sparseProjection.subRowLength[p][i] = subRowPartitionEnd - subRowPartitionBegin;

            // Advance beginning of sub-row partition for next sub-row
            subRowPartitionBegin = subRowPartitionEnd;
        }

        // Add row length to total number of synapses
        s += projection.rowLength(i);
    }
}

template <typename Generator>
SimpleSparseProjection buildFixedProbabilityConnector(unsigned int numPre, unsigned int numPost, float probability, Generator &gen)
{
    // Resize rows
    SimpleSparseProjection projection(numPre);

    // Create RNG to draw probabilities
    std::uniform_real_distribution<> dis(0.0, 1.0);

    // Loop through pre neurons
    for(unsigned int i = 0; i < numPre; i++)
    {
        // Reserve columns based on probability
        // **THINK** might be better to have numPost size vector to generate into
        projection.reserveColumns(i, numPost * probability);

        // Loop through post neurons
        for(unsigned int j = 0; j < numPost; j++)
        {
            // If there should be a connection here, add one to temporary array
            if(dis(gen) < probability)
            {
                projection.addSynapse(i, j);
            }
        }
    }

    return projection;
}