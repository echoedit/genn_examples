#include <cmath>
#include <vector>

#include "modelSpec.h"

#include "../common/exp_curr.h"
#include "../common/lif.h"
#include "../common/vogels_2011.h"


void modelDefinition(NNmodel &model)
{
    initGeNN();
    model.setDT(1.0);
    model.setName("vogels_2011");


    //---------------------------------------------------------------------------
    // Build model
    //---------------------------------------------------------------------------
    // LIF model parameters
    LIF::ParamValues lifParams(
        0.2,    // 0 - C
        20.0,   // 1 - TauM
        -60.0,  // 2 - Vrest
        -60.0,  // 3 - Vreset
        -50.0,  // 4 - Vthresh
        0.2,    // 5 - Ioffset
        5.0);    // 6 - TauRefrac

    // LIF initial conditions
    LIF::VarValues lifInit(
        -55.0,  // 0 - V
        0.0);    // 1 - RefracTime

    // Static synapse parameters
    WeightUpdateModels::StaticPulse::VarValues excitatoryStaticSynapseInit(
        0.03);    // 0 - Wij (nA)

    WeightUpdateModels::StaticPulse::VarValues inhibitoryStaticSynapseInit(
        -0.03);    // 0 - Wij (nA)

    // Additive STDP synapse parameters
    Vogels2011::ParamValues vogels2011AdditiveSTDPParams(
        20.0,   // 0 - Tau
        0.12,   // 1 - rho
        0.005,  // 2 - eta
        -1.0,    // 3 - Wmin
        0.0);    // 4 - Wmax

    Vogels2011::VarValues vogels2011AdditiveSTDPInit(
        0.0);  // 0 - g

    // Exponential current parameters
    ExpCurr::ParamValues excitatoryExpCurrParams(
        5.0);  // 0 - TauSyn (ms)

    ExpCurr::ParamValues inhibitoryExpCurrParams(
        10.0);  // 0 - TauSyn (ms)

    // Create IF_curr neuron
    auto *e = model.addNeuronPopulation<LIF>("E", 2000, lifParams, lifInit);
    model.addNeuronPopulation<LIF>("I", 500, lifParams, lifInit);

    model.addSynapsePopulation<WeightUpdateModels::StaticPulse, ExpCurr>(
        "EE", SynapseMatrixType::SPARSE_GLOBALG, NO_DELAY,
        "E", "E",
        {}, excitatoryStaticSynapseInit,
        excitatoryExpCurrParams, {});
    model.addSynapsePopulation<WeightUpdateModels::StaticPulse, ExpCurr>(
        "EI", SynapseMatrixType::SPARSE_GLOBALG, NO_DELAY,
        "E", "I",
        {}, excitatoryStaticSynapseInit,
        excitatoryExpCurrParams, {});
    model.addSynapsePopulation<WeightUpdateModels::StaticPulse, ExpCurr>(
        "II", SynapseMatrixType::SPARSE_GLOBALG, NO_DELAY,
        "I", "I",
        {}, inhibitoryStaticSynapseInit,
        inhibitoryExpCurrParams, {});
    auto *ie = model.addSynapsePopulation<Vogels2011, ExpCurr>(
        "IE", SynapseMatrixType::SPARSE_INDIVIDUALG, NO_DELAY,
        "I", "E",
        vogels2011AdditiveSTDPParams, vogels2011AdditiveSTDPInit,
        inhibitoryExpCurrParams, {});

    /*auto *ie = model.addSynapsePopulation(
        "IE", NSYNAPSE, SPARSE, INDIVIDUALG, NO_DELAY, MY_EXP_CURR,
        "I", "E",
        staticSynapseInit, NULL,
        NULL, inhibitoryExpCurrParams);*/

    /*model.setSpanTypeToPre("EE");
    model.setSpanTypeToPre("EI");
    model.setSpanTypeToPre("II");
    model.setSpanTypeToPre("IE");*/

    // Use zero-copy for spikes and weights as we want to record them every timestep
    e->setSpikeZeroCopyEnabled(true);
    ie->setWUVarZeroCopyEnabled("g", true);

    model.finalize();
}