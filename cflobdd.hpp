// using namespace CFL_OBDD;

CFLOBDDQuantumCircuit::CFLOBDDQuantumCircuit(unsigned int numQubits, int seed) : QuantumCircuit(numQubits, seed)
{
    // Initialize
        CFLOBDDNodeHandle::InitNoDistinctionTable();
        CFLOBDDNodeHandle::InitAdditionInterleavedTable();
        CFLOBDDNodeHandle::InitReduceCache();
        InitPairProductCache();
        InitTripleProductCache();
        Matrix1234ComplexFloatBoost::Matrix1234Initializer();
        VectorComplexFloatBoost::VectorInitializer();
    //
    unsigned int level = ceil(log2(numQubits)) + 1;
    stateVector = VectorComplexFloatBoost::MkBasisVector(level, 0);
}

CFLOBDDQuantumCircuit::CFLOBDDQuantumCircuit()
{
    // Initialize
        CFLOBDDNodeHandle::InitNoDistinctionTable();
        CFLOBDDNodeHandle::InitAdditionInterleavedTable();
        CFLOBDDNodeHandle::InitReduceCache();
        InitPairProductCache();
        InitTripleProductCache();
        Matrix1234ComplexFloatBoost::Matrix1234Initializer();
        VectorComplexFloatBoost::VectorInitializer();
        numQubits = 0;
    //
}

CFLOBDDQuantumCircuit::~CFLOBDDQuantumCircuit()
{
    DisposeOfTripleProductCache();
	DisposeOfPairProductCache();
	CFLOBDDNodeHandle::DisposeOfReduceCache();
}

void CFLOBDDQuantumCircuit::setNumQubits(unsigned int num)
{
    numQubits = num;
    unsigned int level = ceil(log2(numQubits)) + 1;
    stateVector = VectorComplexFloatBoost::MkBasisVector(level, 0); 
}

CFLOBDD_COMPLEX_BIG ApplyGateF(unsigned int n, unsigned int i, CFLOBDD_COMPLEX_BIG(*f)(unsigned int))
{
    if (n == 1)
    {
        return f(1);
    }
    else {
        int level = ceil(log2(n/2));
        if (i < n/2)
        {
            CFLOBDD_COMPLEX_BIG T = Matrix1234ComplexFloatBoost::MkIdRelationInterleaved(level + 1);
            CFLOBDD_COMPLEX_BIG H = ApplyGateF(n/2, i, f);
            return Matrix1234ComplexFloatBoost::KroneckerProduct2Vocs(H, T);
        }
        else
        {
            CFLOBDD_COMPLEX_BIG T = Matrix1234ComplexFloatBoost::MkIdRelationInterleaved(level + 1);
            return Matrix1234ComplexFloatBoost::KroneckerProduct2Vocs(T, ApplyGateF(n/2, i - n/2, f)); 
        }
    }
}

CFLOBDD_COMPLEX_BIG ApplyGateFWithParam(unsigned int n, unsigned int i, CFLOBDD_COMPLEX_BIG(*f)(unsigned int, double), double theta)
{
    if (n == 1)
    {
        return f(1, theta);
    }
    else {
        int level = ceil(log2(n/2));
        if (i < n/2)
        {
            CFLOBDD_COMPLEX_BIG T = Matrix1234ComplexFloatBoost::MkIdRelationInterleaved(level + 1);
            CFLOBDD_COMPLEX_BIG H = ApplyGateFWithParam(n/2, i, f, theta);
            return Matrix1234ComplexFloatBoost::KroneckerProduct2Vocs(H, T);
        }
        else
        {
            CFLOBDD_COMPLEX_BIG T = Matrix1234ComplexFloatBoost::MkIdRelationInterleaved(level + 1);
            return Matrix1234ComplexFloatBoost::KroneckerProduct2Vocs(T, ApplyGateFWithParam(n/2, i - n/2, f, theta)); 
        }
    }
}

bool checkForInit(unsigned int numQubits)
{
    return numQubits != 0;
}

void CFLOBDDQuantumCircuit::ApplyIdentityGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto H = ApplyGateF(std::pow(2, stateVector.root->level-1), index, Matrix1234ComplexFloatBoost::MkIdRelationInterleaved);
    stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(H, stateVector);
}

void CFLOBDDQuantumCircuit::ApplyHadamardGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto H = ApplyGateF(std::pow(2, stateVector.root->level-1), index, Matrix1234ComplexFloatBoost::MkWalshInterleaved);
    stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(H, stateVector);
    hadamard_count++;
}

void CFLOBDDQuantumCircuit::ApplyNOTGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto X = ApplyGateF(std::pow(2, stateVector.root->level-1), index, Matrix1234ComplexFloatBoost::MkNegationMatrixInterleaved);
    stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(X, stateVector);
}

void CFLOBDDQuantumCircuit::ApplyPauliYGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto Y = ApplyGateF(std::pow(2, stateVector.root->level-1), index, Matrix1234ComplexFloatBoost::MkPauliYMatrixInterleaved);
    stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(Y, stateVector);
}

void CFLOBDDQuantumCircuit::ApplyPauliZGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto Z = ApplyGateF(std::pow(2, stateVector.root->level-1), index, Matrix1234ComplexFloatBoost::MkPauliZMatrixInterleaved);
    stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(Z, stateVector);
}

void CFLOBDDQuantumCircuit::ApplySGate(unsigned int index)
{
    auto S = ApplyGateF(std::pow(2, stateVector.root->level-1), index, Matrix1234ComplexFloatBoost::MkSGateInterleaved);
    stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, stateVector);
}

void CFLOBDDQuantumCircuit::ApplyCNOTGate(long int controller, long int controlled)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(controller != controlled);
    
    if (controller < controlled)
    {
        auto C = Matrix1234ComplexFloatBoost::MkCNOT(stateVector.root->level, std::pow(2, stateVector.root->level - 1), controller, controlled);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else
    {
        auto S = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, controlled, controller);
        auto C = Matrix1234ComplexFloatBoost::MkCNOT(stateVector.root->level, std::pow(2, stateVector.root->level - 1), controlled, controller);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
}

void CFLOBDDQuantumCircuit::ApplySwapGate(long int index1, long int index2)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(index1 != index2);
    
    if (index1 < index2)
    {
        auto C = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, index1, index2);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else
    {
        auto C = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, index2, index1);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
}

void CFLOBDDQuantumCircuit::ApplyiSwapGate(long int index1, long int index2)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(index1 != index2);
    
    if (index1 < index2)
    {
        auto C = Matrix1234ComplexFloatBoost::MkiSwapGate(stateVector.root->level, index1, index2);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else
    {
        auto C = Matrix1234ComplexFloatBoost::MkiSwapGate(stateVector.root->level, index2, index1);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
}

void CFLOBDDQuantumCircuit::ApplyCZGate(long int controller, long int controlled)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(controller != controlled);
    
    if (controller < controlled)
    {
        auto C = Matrix1234ComplexFloatBoost::MkCPGate(stateVector.root->level, controller, controlled, 1.0);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else
    {
        auto S = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, controlled, controller);
        auto C = Matrix1234ComplexFloatBoost::MkCPGate(stateVector.root->level, controlled, controller, 1.0);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
}

void CFLOBDDQuantumCircuit::ApplyCPGate(long int controller, long int controlled, double theta)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(controller != controlled);
    
    if (controller < controlled)
    {
        auto C = Matrix1234ComplexFloatBoost::MkCPGate(stateVector.root->level, controller, controlled, theta);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else
    {
        auto S = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, controlled, controller);
        auto C = Matrix1234ComplexFloatBoost::MkCPGate(stateVector.root->level, controlled, controller, theta);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    } 
}

void CFLOBDDQuantumCircuit::ApplyCSGate(long int controller, long int controlled)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    ApplyCPGate(controller, controlled, 0.5);
}

void CFLOBDDQuantumCircuit::ApplyPhaseShiftGate(unsigned int index, double theta)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto S = ApplyGateFWithParam(std::pow(2, stateVector.root->level-1), index, Matrix1234ComplexFloatBoost::MkPhaseShiftGateInterleaved, theta);
    stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, stateVector); 
}

void CFLOBDDQuantumCircuit::ApplyTGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto S = ApplyGateFWithParam(std::pow(2, stateVector.root->level-1), index, Matrix1234ComplexFloatBoost::MkPhaseShiftGateInterleaved, 0.25);
    stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, stateVector); 
}

void CFLOBDDQuantumCircuit::ApplyCCNOTGate(long int controller1, long int controller2, long int controlled)
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
        auto C = Matrix1234ComplexFloatBoost::MkCCNOT(stateVector.root->level, std::pow(2, stateVector.root->level - 1), controller1, controller2, controlled);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else if (controller1 < controlled && controlled < controller2)
    {
        // a c b   
        auto S = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, controlled, controller2);
        auto C = Matrix1234ComplexFloatBoost::MkCCNOT(stateVector.root->level, std::pow(2, stateVector.root->level - 1), controller1, controlled, controller2);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else if (controller2 < controller1 && controller1 < controlled)
    {
        // b a c
        ApplyCCNOTGate(controller2, controller1, controlled);
    }
    else if (controller2 < controlled && controlled < controller1)
    {
        // b c a
        auto S = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, controlled, controller1);
        auto C = Matrix1234ComplexFloatBoost::MkCCNOT(stateVector.root->level, std::pow(2, stateVector.root->level - 1), controller2, controlled, controller1);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector); 
    }
    else if (controlled < controller1 && controller1 < controller2)
    {
        // c a b
        auto S = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, controlled, controller2);
        // b a c
        auto C = Matrix1234ComplexFloatBoost::MkCCNOT(stateVector.root->level, std::pow(2, stateVector.root->level - 1), controlled, controller1, controller2);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else if (controlled < controller2 && controller2 < controller1)
    {
        // c b a
        auto S = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, controlled, controller1);
        // a b c
        auto C = Matrix1234ComplexFloatBoost::MkCCNOT(stateVector.root->level, std::pow(2, stateVector.root->level - 1), controlled, controller2, controller1);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
}

void CFLOBDDQuantumCircuit::ApplyCSwapGate(long int controller, long int index1, long int index2)
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
        auto C = Matrix1234ComplexFloatBoost::MkCSwapGate(stateVector.root->level, controller, index1, index2);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else if (controller < index2 && index2 < index1)
    {
        // a c b   
        auto C = Matrix1234ComplexFloatBoost::MkCSwapGate(stateVector.root->level, controller, index2, index1);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else if (index1 < controller && controller < index2)
    {
        // b a c
        auto S = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, index1, controller);
        auto C = Matrix1234ComplexFloatBoost::MkCSwapGate(stateVector.root->level, index1, controller, index2);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else if (index1 < index2 && index2 < controller)
    {
        // b c a
        auto S = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, index1, controller);
        auto C = Matrix1234ComplexFloatBoost::MkCSwapGate(stateVector.root->level, index1, index2, controller);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else if (index2 < controller && controller < index1)
    {
        // c a b
        auto S = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, index2, controller);
        // b a c
        auto C = Matrix1234ComplexFloatBoost::MkCSwapGate(stateVector.root->level, index2, controller, index1);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else if (index2 < index1 && index1 < controller)
    {
        // c b a
        auto S = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, index2, controller);
        // a b c
        auto C = Matrix1234ComplexFloatBoost::MkCSwapGate(stateVector.root->level, index2, index1, controller);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
}

void CFLOBDDQuantumCircuit::ApplyGlobalPhase(double phase)
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

long double CFLOBDDQuantumCircuit::GetProbability(std::map<unsigned int, int>& qubit_vals)
{
    // stateVector.print(std::cout);
    auto tmp = VectorComplexFloatBoost::VectorWithAmplitude(stateVector);
    std::string s(std::pow(2, tmp.root->level-1), 'X');
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
    auto restricted = Matrix1234ComplexFloatBoost::MkRestrictMatrix(tmp.root->level, s);
    tmp = tmp * restricted;
    tmp.CountPaths();
    return VectorComplexFloatBoost::getNonZeroProbability(tmp);
}

std::string CFLOBDDQuantumCircuit::Measure() 
{
    auto tmp = VectorComplexFloatBoost::VectorWithAmplitude(stateVector);
    tmp.CountPaths();
    return VectorComplexFloatBoost::Sampling(tmp, true).substr(0, numQubits); 
}

unsigned long long int CFLOBDDQuantumCircuit::GetPathCount(long double prob)
{
    auto tmp = VectorComplexFloatBoost::VectorWithAmplitude(stateVector);
    tmp.CountPaths();
    return VectorComplexFloatBoost::GetPathCount(tmp, prob);  
}
