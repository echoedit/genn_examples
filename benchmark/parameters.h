#pragma once

//------------------------------------------------------------------------
// Parameters
//------------------------------------------------------------------------
#define SYNAPSE_MATRIX_SPARSE
#define SYNAPSE_MATRIX_INDIVIDUAL
//#define VALIDATE

namespace Parameters
{
    const unsigned int numPre = 30000;
    const unsigned int numPost = 30000;

    const float connectionProbability = 0.1f;
    const float weight = 0.001f;

#ifdef SYNAPSE_MATRIX_SPARSE
#ifdef SYNAPSE_MATRIX_INDIVIDUAL
    SynapseMatrixType synapseMatrixType = SynapseMatrixType::SPARSE_INDIVIDUALG;
#else
    SynapseMatrixType synapseMatrixType = SynapseMatrixType::SPARSE_GLOBALG;
#endif
#else   // SYNAPSE_MATRIX_SPARSE
    SynapseMatrixType synapseMatrixType = SynapseMatrixType::DENSE_INDIVIDUALG;
#endif
}