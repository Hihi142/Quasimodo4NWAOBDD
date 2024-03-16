
// *******************
// Weighted BDD
// *******************

#include "cflobdd/CFLOBDD/wvector_complex_fb_mul.h"
#include "cflobdd/CFLOBDD/weighted_cross_product.h"
#include "cflobdd/CFLOBDD/weighted_cross_product_bdd.h"

WeightedBDDQuantumCircuit::WeightedBDDQuantumCircuit(unsigned int numQubits, int seed) : QuantumCircuit(numQubits, seed)
{
    // Initialize
    WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::InitNoDistinctionTable();
	WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::InitNoDistinctionTable_Ann();
	WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::InitIdentityNodeTable();	
	WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::InitReduceCache();
	WeightedMatrix1234ComplexFloatBoostMul::Matrix1234Initializer();
	WeightedVectorComplexFloatBoostMul::VectorInitializer();
	InitWeightedPairProductCache<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>();

    WeightedBDDNodeHandle<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::InitLeafNodes();
	InitWeightedBDDPairProductCache<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>();
    //
    std::string s(2 * numQubits, '0');
    stateVector = WeightedVectorComplexFloatBoostMul::MkBasisVector(2 * numQubits, s, 0);
}

WeightedBDDQuantumCircuit::WeightedBDDQuantumCircuit()
{
    // Initialize
    WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::InitNoDistinctionTable();
	WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::InitNoDistinctionTable_Ann();
	WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::InitIdentityNodeTable();	
	WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::InitReduceCache();
	WeightedMatrix1234ComplexFloatBoostMul::Matrix1234Initializer();
	WeightedVectorComplexFloatBoostMul::VectorInitializer();
	InitWeightedPairProductCache<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>();

    WeightedBDDNodeHandle<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::InitLeafNodes();
	InitWeightedBDDPairProductCache<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>();
    numQubits = 0;
    //
}

WeightedBDDQuantumCircuit::~WeightedBDDQuantumCircuit()
{
    DisposeOfWeightedPairProductCache<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>();
    DisposeOfWeightedBDDPairProductCache<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>();
}

void WeightedBDDQuantumCircuit::setNumQubits(unsigned int num)
{
    numQubits = num;
    stateVector = WeightedVectorComplexFloatBoostMul::MkBasisVector(numQubits, 0, 0);
}

WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL ApplyGateF(unsigned int n, unsigned int i, WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(*f)(unsigned int, int))
{
    WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL H = f(2, 0);
    if (i == 0)
    {
        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL I = WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved(2 * (n-1), 0);
        return WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(H, I); 
    }
    else if (i == n - 1)
    {
        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL I = WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved(2 * (n-1), 0);
        return WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(I, H);
    }
    else {
        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL I1 = WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved(2 * i, 0);
        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL I2 = WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved(2 * (n - i - 1), 0);
        auto T = WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(H, I2);
        return WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(I1, T);
    }
}

WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL ApplyGateFWithParam(unsigned int n, unsigned int i, WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(*f)(unsigned int, double, int), double theta)
{
    WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL H = f(2, theta, 0);
    if (i == 0)
    {
        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL I = WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved(2 * (n-1), 0);
        return WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(H, I); 
    }
    else if (i == n - 1)
    {
        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL I = WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved(2 * (n-1), 0);
        return WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(I, H);
    }
    else {
        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL I1 = WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved(2 * i, 0);
        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL I2 = WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved(2 * (n - i - 1), 0);
        auto T = WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(H, I2);
        return WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(I1, T);
    }
}


void WeightedBDDQuantumCircuit::ApplyIdentityGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto H = ApplyGateF(numQubits, index, WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved);
    stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(H, stateVector);
}

void WeightedBDDQuantumCircuit::ApplyHadamardGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto H = ApplyGateF(numQubits, index, WeightedMatrix1234ComplexFloatBoostMul::MkWalshInterleaved);
    // H.print(std::cout);
    stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(H, stateVector);
    // stateVector.print(std::cout);
    // std::cout << std::endl;
    hadamard_count++;
}

void WeightedBDDQuantumCircuit::ApplyNOTGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    // stateVector.print(std::cout);
    auto X = ApplyGateF(numQubits, index, WeightedMatrix1234ComplexFloatBoostMul::MkNegationMatrixInterleaved);
    stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(X, stateVector);
    // stateVector.print(std::cout);
}

void WeightedBDDQuantumCircuit::ApplyPauliYGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto Y = ApplyGateF(numQubits, index, WeightedMatrix1234ComplexFloatBoostMul::MkPauliYGate);
    stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(Y, stateVector);
}

void WeightedBDDQuantumCircuit::ApplyPauliZGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto Z = ApplyGateF(numQubits, index, WeightedMatrix1234ComplexFloatBoostMul::MkPauliZGate);
    stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(Z, stateVector);
}

void WeightedBDDQuantumCircuit::ApplySGate(unsigned int index)
{
    auto S = ApplyGateF(numQubits, index, WeightedMatrix1234ComplexFloatBoostMul::MkSGate);
    stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, stateVector);
}

void WeightedBDDQuantumCircuit::ApplyCNOTGate(long int controller, long int controlled)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(controller != controlled);
    
    if (controller < controlled)
    {
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCNOT(2 * numQubits, numQubits, controller, controlled, 0);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
        // stateVector.print(std::cout);
    }
    else
    {
        auto S = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, controlled, controller, 0);
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCNOT(2 * numQubits, numQubits, controlled, controller, 0);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, S);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, C);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
}

void WeightedBDDQuantumCircuit::ApplySwapGate(long int index1, long int index2)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(index1 != index2);
    
    if (index1 < index2)
    {
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, index1, index2, 0);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else
    {
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, index2, index1, 0);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
}

void WeightedBDDQuantumCircuit::ApplyiSwapGate(long int index1, long int index2)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(index1 != index2);
    
    if (index1 < index2)
    {
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkiSwapGate(2 * numQubits, index1, index2, 0);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else
    {
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkiSwapGate(2 * numQubits, index2, index1, 0);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
}

void WeightedBDDQuantumCircuit::ApplyCZGate(long int controller, long int controlled)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(controller != controlled);
    
    if (controller < controlled)
    {
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCPGate(2 * numQubits, controller, controlled, 1.0, 0);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else
    {
        auto S = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, controlled, controller, 0);
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCPGate(2 * numQubits, controlled, controller, 1.0, 0);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, S);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, C);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
}

void WeightedBDDQuantumCircuit::ApplyCPGate(long int controller, long int controlled, double theta)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(controller != controlled);
    
    if (controller < controlled)
    {
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCPGate(2 * numQubits, controller, controlled, theta, 0);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else
    {
        auto S = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, controlled, controller, 0);
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCPGate(2 * numQubits, controlled, controller, theta, 0);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, S);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, C);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    } 
}

void WeightedBDDQuantumCircuit::ApplyCSGate(long int controller, long int controlled)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    ApplyCPGate(controller, controlled, 0.5);
}

void WeightedBDDQuantumCircuit::ApplyPhaseShiftGate(unsigned int index, double theta)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto S = ApplyGateFWithParam(numQubits, index, WeightedMatrix1234ComplexFloatBoostMul::MkPhaseShiftGate, theta);
    stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, stateVector); 
}

void WeightedBDDQuantumCircuit::ApplyTGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto S = ApplyGateFWithParam(numQubits, index, WeightedMatrix1234ComplexFloatBoostMul::MkPhaseShiftGate, 0.25);
    stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, stateVector); 
}

void WeightedBDDQuantumCircuit::ApplyCCNOTGate(long int controller1, long int controller2, long int controlled)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(controller1 != controlled);
    assert(controller2 != controlled);
    assert(controller1 != controller2);
    if (controller1 < controller2 && controller2 < controlled)
    {
        // a b c
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCCNOT(2 * numQubits, controller1, controller2, controlled, 0);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else if (controller1 < controlled && controlled < controller2)
    {
        // a c b   
        auto S = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, controlled, controller2, 0);
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCCNOT(2 * numQubits, controller1, controlled, controller2, 0);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, S);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, C);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else if (controller2 < controller1 && controller1 < controlled)
    {
        // b a c
        ApplyCCNOTGate(controller2, controller1, controlled);
    }
    else if (controller2 < controlled && controlled < controller1)
    {
        // b c a
        auto S = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, controlled, controller1, 0);
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCCNOT(2 * numQubits, controller2, controlled, controller1, 0);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, S);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, C);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector); 
    }
    else if (controlled < controller1 && controller1 < controller2)
    {
        // c a b
        auto S = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, controlled, controller2, 0);
        // b a c
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCCNOT(2 * numQubits, controlled, controller1, controller2, 0);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, S);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, C);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else if (controlled < controller2 && controller2 < controller1)
    {
        // c b a
        auto S = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, controlled, controller1, 0);
        // a b c
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCCNOT(2 * numQubits, controlled, controller2, controller1, 0);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, S);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, C);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
}

void WeightedBDDQuantumCircuit::ApplyCSwapGate(long int controller, long int index1, long int index2)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(controller != index1);
    assert(controller != index2);
    assert(index1 != index2);
    
    if (controller < index1 && index1 < index2)
    {
        // a b c
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCSwapGate(2 * numQubits, controller, index1, index2, 0);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else if (controller < index2 && index2 < index1)
    {
        // a c b   
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCSwapGate(2 * numQubits, controller, index2, index1, 0);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else if (index1 < controller && controller < index2)
    {
        // b a c
        auto S = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, index1, controller, 0);
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCSwapGate(2 * numQubits, index1, controller, index2, 0);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, S);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, C);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else if (index1 < index2 && index2 < controller)
    {
        // b c a
        auto S = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, index1, controller, 0);
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCSwapGate(2 * numQubits, index1, index2, controller, 0);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, S);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, C);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else if (index2 < controller && controller < index1)
    {
        // c a b
        auto S = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, index2, controller, 0);
        // b a c
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCSwapGate(2 * numQubits, index2, controller, index1, 0);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, S);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, C);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else if (index2 < index1 && index1 < controller)
    {
        // c b a
        auto S = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, index2, controller, 0);
        // a b c
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCSwapGate(2 * numQubits, index2, index1, controller, 0);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, S);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, C);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
}

void WeightedBDDQuantumCircuit::ApplyGlobalPhase(double phase)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto cos_v = boost::math::cos_pi(phase);
    auto sin_v = boost::math::sin_pi(phase);
    BIG_COMPLEX_FLOAT phase_complex(cos_v, sin_v);
    stateVector = phase_complex * stateVector;
}

long double WeightedBDDQuantumCircuit::GetProbability(std::map<unsigned int, int>& qubit_vals)
{
    auto tmp = stateVector;
    std::string s(numQubits, 'X');
    for (unsigned int i = 0; i < numQubits; i++)
    {
        if (qubit_vals.find(i) != qubit_vals.end())
        {
            if (qubit_vals[i] == 0)
                s[i] = '0';
            else if (qubit_vals[i] == 1)
                s[i] = '1';   
        }
    }
    auto restricted = WeightedMatrix1234ComplexFloatBoostMul::MkRestrictMatrix(2 * numQubits, s, 0);
    tmp = tmp * restricted;
    tmp.ComputeWeightOfPathsAsAmpsToExits();
    return WeightedVectorComplexFloatBoostMul::getNonZeroProbability(tmp);
}

std::string WeightedBDDQuantumCircuit::Measure() 
{
    auto tmp = stateVector;
    // tmp.print(std::cout);
    tmp.ComputeWeightOfPathsAsAmpsToExits();
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    return WeightedVectorComplexFloatBoostMul::Sampling(tmp, true, mt, dis).substr(0, numQubits); 
}

unsigned long long int WeightedBDDQuantumCircuit::GetPathCount(long double prob)
{
    std::cout << "Error! Operation not supported in WBDDs" << std::endl;
    abort();
} 

