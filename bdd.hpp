
/// ****** BDDQuantumCircuit *******

#include <mpfr.h>
#define RND_TYPE MPFR_RNDF

BDDQuantumCircuit::BDDQuantumCircuit(unsigned int numQubits, int seed) : QuantumCircuit(numQubits, seed)
{
    mgr = new Cudd(0,0);

    if (numQubits > 512) // Based on experiments
    {
        mpfr_set_default_prec(300);
        CUDD_VALUE_TYPE epsilon;
        mpfr_init_set_si(epsilon.real, -1 * (200) , RND_TYPE);
        mpfr_exp10(epsilon.real, epsilon.real, RND_TYPE);
        mpfr_init_set_si(epsilon.imag, 0, RND_TYPE);
        mgr->SetEpsilon(epsilon);
    }

    if (numQubits > 2048) // Based on experiments
    {
        mpfr_set_default_prec(500);
        CUDD_VALUE_TYPE epsilon;
        mpfr_init_set_si(epsilon.real, -1 * (300) , RND_TYPE);
        mpfr_exp10(epsilon.real, epsilon.real, RND_TYPE);
        mpfr_init_set_si(epsilon.imag, 0, RND_TYPE);
        mgr->SetEpsilon(epsilon);
    }

    for (unsigned int i = 0; i < numQubits; i++)
    {
        x_vars.push_back(mgr->addVar(2*i));
        y_vars.push_back(mgr->addVar(2*i + 1));
    }

    stateVector = mgr->addOne();
    // e_{0..0}
    for (unsigned int i = 0; i < numQubits; i++)
    {
        stateVector *= ~x_vars[i];
    }
}

BDDQuantumCircuit::BDDQuantumCircuit()
{
    mgr = new Cudd(0,0);
    numQubits = 0;
}

BDDQuantumCircuit::~BDDQuantumCircuit()
{
    // delete mgr;
}

void BDDQuantumCircuit::setNumQubits(unsigned int n)
{
    numQubits = n;
    x_vars.clear();
    y_vars.clear();
    for (unsigned int i = 0; i < numQubits; i++)
    {
        x_vars.push_back(mgr->addVar(2*i));
        y_vars.push_back(mgr->addVar(2*i + 1));
    }

    stateVector = mgr->addOne();
    // e_{0..0}
    for (unsigned int i = 0; i < numQubits; i++)
    {
        stateVector *= ~x_vars[i];
    }
}

void BDDQuantumCircuit::ApplyIdentityGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    ADD IDGate = ~y_vars[index] * ~x_vars[index] + y_vars[index] * x_vars[index];
    std::vector<ADD> tmp_x, tmp_y; 
    tmp_x.push_back(x_vars[index]); tmp_y.push_back(y_vars[index]);
    stateVector = IDGate.MatrixMultiply(stateVector, tmp_x);
    stateVector = stateVector.SwapVariables(tmp_y, tmp_x);
}

void BDDQuantumCircuit::ApplyHadamardGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    
    CUDD_VALUE_TYPE H;
    mpfr_init_set_d(H.real, 1, RND_TYPE);
    mpfr_init_set_d(H.imag, 0, RND_TYPE);
    mpfr_div_d(H.real, H.real, sqrt(2), RND_TYPE);
    ADD H_val = mgr->constant(H);
    mpfr_clear(H.real); mpfr_clear(H.imag);
    
    ADD HGate = (~y_vars[index] + y_vars[index] * (~x_vars[index] - x_vars[index])) * H_val;
    // HGate.print(2*numQubits, 2);
    std::vector<ADD> tmp_x, tmp_y; 
    tmp_x.push_back(x_vars[index]); tmp_y.push_back(y_vars[index]);
    stateVector = HGate.MatrixMultiply(stateVector, tmp_x);
    stateVector = stateVector.SwapVariables(tmp_y, tmp_x);
    // if (index == numQubits - 1)
    //     stateVector.print(2 * numQubits, 2);
    hadamard_count++;
}

void BDDQuantumCircuit::ApplyNOTGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    ADD XGate = ~y_vars[index] * x_vars[index] + y_vars[index] * ~x_vars[index];
    std::vector<ADD> tmp_x, tmp_y; 
    tmp_x.push_back(x_vars[index]); tmp_y.push_back(y_vars[index]);
    // std::cout << "index: " << index << std::endl;
    // stateVector.print(2 * numQubits, 2);
    stateVector = XGate.MatrixMultiply(stateVector, tmp_x);
    stateVector = stateVector.SwapVariables(tmp_y, tmp_x);
    // std::cout << "after" << std::endl;
    // stateVector.print(2 * numQubits, 2);
}

void BDDQuantumCircuit::ApplyPauliYGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    CUDD_VALUE_TYPE I;
    mpfr_init_set_d(I.real, 0, RND_TYPE);
    mpfr_init_set_d(I.imag, 1, RND_TYPE);
    ADD I_val = mgr->constant(I);
    mpfr_clear(I.real); mpfr_clear(I.imag);

    ADD Y =  (y_vars[index] * x_vars[index] - ~y_vars[index] * x_vars[index]) * I_val;
    std::vector<ADD> tmp_x, tmp_y; 
    tmp_x.push_back(x_vars[index]); tmp_y.push_back(y_vars[index]);
    stateVector = Y.MatrixMultiply(stateVector, tmp_x);
    stateVector = stateVector.SwapVariables(tmp_y, tmp_x);
}

void BDDQuantumCircuit::ApplyPauliZGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    ADD ZGate = ~y_vars[index] * ~x_vars[index] - y_vars[index] * x_vars[index];
    std::vector<ADD> tmp_x, tmp_y; 
    tmp_x.push_back(x_vars[index]); tmp_y.push_back(y_vars[index]);
    stateVector = ZGate.MatrixMultiply(stateVector, tmp_x);
    stateVector = stateVector.SwapVariables(tmp_y, tmp_x);
}

void BDDQuantumCircuit::ApplySGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    CUDD_VALUE_TYPE I;
    mpfr_init_set_d(I.real, 0, RND_TYPE);
    mpfr_init_set_d(I.imag, 1, RND_TYPE);
    ADD I_val = mgr->constant(I);
    mpfr_clear(I.real); mpfr_clear(I.imag);

    ADD SGate = ~y_vars[index] * ~x_vars[index] + y_vars[index] * x_vars[index] * I_val;
    std::vector<ADD> tmp_x, tmp_y; 
    tmp_x.push_back(x_vars[index]); tmp_y.push_back(y_vars[index]);
    stateVector = SGate.MatrixMultiply(stateVector, tmp_x);
    stateVector = stateVector.SwapVariables(tmp_y, tmp_x);
}

void BDDQuantumCircuit::ApplyCNOTGate(long int controller, long int controlled)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    ADD CNOTGate = ~y_vars[controller] * ~y_vars[controlled] * ~x_vars[controller] * ~x_vars[controlled]
                 + ~y_vars[controller] * y_vars[controlled] * ~x_vars[controller] * x_vars[controlled]
                 + y_vars[controller] * ~y_vars[controlled] * x_vars[controller] * x_vars[controlled]
                 + y_vars[controller] * y_vars[controlled] * x_vars[controller] * ~x_vars[controlled];
    std::vector<ADD> tmp_x, tmp_y; 
    tmp_x.push_back(x_vars[controller]); tmp_x.push_back(x_vars[controlled]);
    tmp_y.push_back(y_vars[controller]); tmp_y.push_back(y_vars[controlled]);
    stateVector = CNOTGate.MatrixMultiply(stateVector, tmp_x);
    stateVector = stateVector.SwapVariables(tmp_y, tmp_x);
}

void BDDQuantumCircuit::ApplyGlobalPhase(double phase)
{
   if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    CUDD_VALUE_TYPE phase_complex;
    mpfr_init(phase_complex.real); mpfr_init(phase_complex.imag);
    mpfr_const_pi(phase_complex.real, RND_TYPE);
    mpfr_const_pi(phase_complex.imag, RND_TYPE);
    mpfr_mul_d(phase_complex.real, phase_complex.real, phase, RND_TYPE);
    mpfr_mul_d(phase_complex.imag, phase_complex.imag, phase, RND_TYPE);
    mpfr_cos(phase_complex.real, phase_complex.real, RND_TYPE);
    mpfr_sin(phase_complex.imag, phase_complex.imag, RND_TYPE);
    ADD phase_add = mgr->constant(phase_complex);
    mpfr_clear(phase_complex.real); mpfr_clear(phase_complex.imag);
    stateVector = phase_add * stateVector; 
}

void BDDQuantumCircuit::ApplySwapGate(long int index1, long int index2)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    ADD SwapGate = ~y_vars[index1] * ~y_vars[index2] * ~x_vars[index1] * ~x_vars[index2]
                 + ~y_vars[index1] * y_vars[index2] * x_vars[index1] * ~x_vars[index2]
                 + y_vars[index1] * ~y_vars[index2] * ~x_vars[index1] * x_vars[index2]
                 + y_vars[index1] * y_vars[index2] * x_vars[index1] * x_vars[index2];
    std::vector<ADD> tmp_x, tmp_y; 
    tmp_x.push_back(x_vars[index1]); tmp_x.push_back(x_vars[index2]);
    tmp_y.push_back(y_vars[index1]); tmp_y.push_back(y_vars[index2]);
    stateVector = SwapGate.MatrixMultiply(stateVector, tmp_x);
    stateVector = stateVector.SwapVariables(tmp_y, tmp_x);
}

void BDDQuantumCircuit::ApplyiSwapGate(long int index1, long int index2)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }

    CUDD_VALUE_TYPE I;
    mpfr_init_set_d(I.real, 0, RND_TYPE);
    mpfr_init_set_d(I.imag, 1, RND_TYPE);
    ADD I_val = mgr->constant(I);
    mpfr_clear(I.real); mpfr_clear(I.imag);

    ADD iSwapGate = ~y_vars[index1] * ~y_vars[index2] * ~x_vars[index1] * ~x_vars[index2]
                 + ~y_vars[index1] * y_vars[index2] * x_vars[index1] * ~x_vars[index2] * I_val
                 + y_vars[index1] * ~y_vars[index2] * ~x_vars[index1] * x_vars[index2] * I_val
                 + y_vars[index1] * y_vars[index2] * x_vars[index1] * x_vars[index2];
    std::vector<ADD> tmp_x, tmp_y; 
    tmp_x.push_back(x_vars[index1]); tmp_x.push_back(x_vars[index2]);
    tmp_y.push_back(y_vars[index1]); tmp_y.push_back(y_vars[index2]);
    stateVector = iSwapGate.MatrixMultiply(stateVector, tmp_x);
    stateVector = stateVector.SwapVariables(tmp_y, tmp_x);
}

void BDDQuantumCircuit::ApplyCZGate(long int controller, long int controlled)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    ADD CZGate = ~y_vars[controller] * ~y_vars[controlled] * ~x_vars[controller] * ~x_vars[controlled]
                 + ~y_vars[controller] * y_vars[controlled] * ~x_vars[controller] * x_vars[controlled]
                 + y_vars[controller] * ~y_vars[controlled] * x_vars[controller] * ~x_vars[controlled]
                 - y_vars[controller] * y_vars[controlled] * x_vars[controller] * x_vars[controlled];
    std::vector<ADD> tmp_x, tmp_y; 
    tmp_x.push_back(x_vars[controller]); tmp_x.push_back(x_vars[controlled]);
    tmp_y.push_back(y_vars[controller]); tmp_y.push_back(y_vars[controlled]);
    stateVector = CZGate.MatrixMultiply(stateVector, tmp_x);
    stateVector = stateVector.SwapVariables(tmp_y, tmp_x);
}

void BDDQuantumCircuit::ApplyCPGate(long int controller, long int controlled, double theta)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }

    CUDD_VALUE_TYPE phase_complex;
    mpfr_init(phase_complex.real); mpfr_init(phase_complex.imag);
    mpfr_const_pi(phase_complex.real, RND_TYPE);
    mpfr_const_pi(phase_complex.imag, RND_TYPE);
    mpfr_mul_d(phase_complex.real, phase_complex.real, theta, RND_TYPE);
    mpfr_mul_d(phase_complex.imag, phase_complex.imag, theta, RND_TYPE);
    mpfr_cos(phase_complex.real, phase_complex.real, RND_TYPE);
    mpfr_sin(phase_complex.imag, phase_complex.imag, RND_TYPE);
    ADD phase_add = mgr->constant(phase_complex);
    mpfr_clear(phase_complex.real); mpfr_clear(phase_complex.imag);


    ADD CPGate = ~y_vars[controller] * ~y_vars[controlled] * ~x_vars[controller] * ~x_vars[controlled]
                 + ~y_vars[controller] * y_vars[controlled] * ~x_vars[controller] * x_vars[controlled]
                 + y_vars[controller] * ~y_vars[controlled] * x_vars[controller] * x_vars[controlled]
                 + y_vars[controller] * y_vars[controlled] * x_vars[controller] * ~x_vars[controlled] * phase_add;
    std::vector<ADD> tmp_x, tmp_y; 
    tmp_x.push_back(x_vars[controller]); tmp_x.push_back(x_vars[controlled]);
    tmp_y.push_back(y_vars[controller]); tmp_y.push_back(y_vars[controlled]);
    stateVector = CPGate.MatrixMultiply(stateVector, tmp_x);
    stateVector = stateVector.SwapVariables(tmp_y, tmp_x);
}

void BDDQuantumCircuit::ApplyPhaseShiftGate(unsigned int index, double theta)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    CUDD_VALUE_TYPE phase_complex;
    mpfr_init(phase_complex.real); mpfr_init(phase_complex.imag);
    mpfr_const_pi(phase_complex.real, RND_TYPE);
    mpfr_const_pi(phase_complex.imag, RND_TYPE);
    mpfr_mul_d(phase_complex.real, phase_complex.real, theta, RND_TYPE);
    mpfr_mul_d(phase_complex.imag, phase_complex.imag, theta, RND_TYPE);
    mpfr_cos(phase_complex.real, phase_complex.real, RND_TYPE);
    mpfr_sin(phase_complex.imag, phase_complex.imag, RND_TYPE);
    ADD phase_add = mgr->constant(phase_complex);
    mpfr_clear(phase_complex.real); mpfr_clear(phase_complex.imag);

    ADD PhaseShiftGate =  ~y_vars[index] * ~x_vars[index] + y_vars[index] * x_vars[index] * phase_add;
    std::vector<ADD> tmp_x, tmp_y; 
    tmp_x.push_back(x_vars[index]); tmp_y.push_back(y_vars[index]);
    stateVector = PhaseShiftGate.MatrixMultiply(stateVector, tmp_x);
    stateVector = stateVector.SwapVariables(tmp_y, tmp_x);
}

void BDDQuantumCircuit::ApplyTGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    CUDD_VALUE_TYPE phase_complex;
    mpfr_init(phase_complex.real); mpfr_init(phase_complex.imag);
    mpfr_const_pi(phase_complex.real, RND_TYPE);
    mpfr_const_pi(phase_complex.imag, RND_TYPE);
    mpfr_mul_d(phase_complex.real, phase_complex.real, 0.25, RND_TYPE);
    mpfr_mul_d(phase_complex.imag, phase_complex.imag, 0.25, RND_TYPE);
    mpfr_cos(phase_complex.real, phase_complex.real, RND_TYPE);
    mpfr_sin(phase_complex.imag, phase_complex.imag, RND_TYPE);
    ADD phase_add = mgr->constant(phase_complex);
    mpfr_clear(phase_complex.real); mpfr_clear(phase_complex.imag);

    ADD TGate =  ~y_vars[index] * ~x_vars[index] + y_vars[index] * x_vars[index] * phase_add;
    std::vector<ADD> tmp_x, tmp_y; 
    tmp_x.push_back(x_vars[index]); tmp_y.push_back(y_vars[index]); 
    stateVector = TGate.MatrixMultiply(stateVector, tmp_x);
    stateVector = stateVector.SwapVariables(tmp_y, tmp_x);
}

void BDDQuantumCircuit::ApplyCSGate(long int controller, long int controlled)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }

    CUDD_VALUE_TYPE phase_complex;
    mpfr_init(phase_complex.real); mpfr_init(phase_complex.imag);
    mpfr_const_pi(phase_complex.real, RND_TYPE);
    mpfr_const_pi(phase_complex.imag, RND_TYPE);
    mpfr_mul_d(phase_complex.real, phase_complex.real, 0.5, RND_TYPE);
    mpfr_mul_d(phase_complex.imag, phase_complex.imag, 0.5, RND_TYPE);
    mpfr_cos(phase_complex.real, phase_complex.real, RND_TYPE);
    mpfr_sin(phase_complex.imag, phase_complex.imag, RND_TYPE);
    ADD phase_add = mgr->constant(phase_complex);
    mpfr_clear(phase_complex.real); mpfr_clear(phase_complex.imag);

    ADD CSGate = ~y_vars[controller] * ~y_vars[controlled] * ~x_vars[controller] * ~x_vars[controlled]
                 + ~y_vars[controller] * y_vars[controlled] * ~x_vars[controller] * x_vars[controlled]
                 + y_vars[controller] * ~y_vars[controlled] * x_vars[controller] * ~x_vars[controlled]
                 + y_vars[controller] * y_vars[controlled] * x_vars[controller] * x_vars[controlled] * phase_add;
    std::vector<ADD> tmp_x, tmp_y; 
    tmp_x.push_back(x_vars[controller]); tmp_x.push_back(x_vars[controlled]);
    tmp_y.push_back(y_vars[controller]); tmp_y.push_back(y_vars[controlled]);
    stateVector = CSGate.MatrixMultiply(stateVector, tmp_x);
    stateVector = stateVector.SwapVariables(tmp_y, tmp_x);
}

void BDDQuantumCircuit::ApplyCCNOTGate(long int controller, long int index1, long int index2)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    
    ADD CCNOTGate = ~y_vars[controller] * ~y_vars[index1] * ~y_vars[index2] * ~x_vars[controller] * ~x_vars[index1] * ~x_vars[index2]
                  + ~y_vars[controller] * ~y_vars[index1] * y_vars[index2] * ~x_vars[controller] * ~x_vars[index1] * x_vars[index2]
                  + ~y_vars[controller] * y_vars[index1] * ~y_vars[index2] * ~x_vars[controller] * x_vars[index1] * ~x_vars[index2]
                  + ~y_vars[controller] * y_vars[index1] * y_vars[index2] * ~x_vars[controller] * x_vars[index1] * x_vars[index2]
                  + y_vars[controller] * ~y_vars[index1] * ~y_vars[index2] * x_vars[controller] * ~x_vars[index1] * ~x_vars[index2]
                  + y_vars[controller] * ~y_vars[index1] * y_vars[index2] * x_vars[controller] * ~x_vars[index1] * x_vars[index2]
                  + y_vars[controller] * y_vars[index1] * ~y_vars[index2] * x_vars[controller] * x_vars[index1] * x_vars[index2]
                  + y_vars[controller] * y_vars[index1] * y_vars[index2] * x_vars[controller] * x_vars[index1] * ~x_vars[index2];
    std::vector<ADD> tmp_x, tmp_y; 
    tmp_x.push_back(x_vars[controller]); tmp_x.push_back(x_vars[index1]); tmp_x.push_back(x_vars[index2]);
    tmp_y.push_back(y_vars[controller]); tmp_y.push_back(y_vars[index1]); tmp_y.push_back(y_vars[index2]);
    stateVector = CCNOTGate.MatrixMultiply(stateVector, tmp_x);
    stateVector = stateVector.SwapVariables(tmp_y, tmp_x);
}

void BDDQuantumCircuit::ApplyCSwapGate(long int controller, long int index1, long int index2)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    ADD CSwapGate = ~y_vars[controller] * ~y_vars[index1] * ~y_vars[index2] * ~x_vars[controller] * ~x_vars[index1] * ~x_vars[index2]
                  + ~y_vars[controller] * ~y_vars[index1] * y_vars[index2] * ~x_vars[controller] * ~x_vars[index1] * x_vars[index2]
                  + ~y_vars[controller] * y_vars[index1] * ~y_vars[index2] * ~x_vars[controller] * x_vars[index1] * ~x_vars[index2]
                  + ~y_vars[controller] * y_vars[index1] * y_vars[index2] * ~x_vars[controller] * x_vars[index1] * x_vars[index2]
                  + y_vars[controller] * ~y_vars[index1] * ~y_vars[index2] * x_vars[controller] * ~x_vars[index1] * ~x_vars[index2]
                  + y_vars[controller] * ~y_vars[index1] * y_vars[index2] * x_vars[controller] * x_vars[index1] * ~x_vars[index2]
                  + y_vars[controller] * y_vars[index1] * ~y_vars[index2] * x_vars[controller] * ~x_vars[index1] * x_vars[index2]
                  + y_vars[controller] * y_vars[index1] * y_vars[index2] * x_vars[controller] * x_vars[index1] * x_vars[index2];
    std::vector<ADD> tmp_x, tmp_y; 
    tmp_x.push_back(x_vars[controller]); tmp_x.push_back(x_vars[index1]); tmp_x.push_back(x_vars[index2]);
    tmp_y.push_back(y_vars[controller]); tmp_y.push_back(y_vars[index1]); tmp_y.push_back(y_vars[index2]);
    stateVector = CSwapGate.MatrixMultiply(stateVector, tmp_x);
    stateVector = stateVector.SwapVariables(tmp_y, tmp_x);
}

long double BDDQuantumCircuit::GetProbability(std::map<unsigned int, int>& qubit_vals)
{
    // stateVector.print(2*numQubits, 2);
    ADD tmp = stateVector.SquareTerminalValues();
    ADD s_add = mgr->addOne();

    for (auto it = qubit_vals.begin(); it != qubit_vals.end(); it++)
    {
        if (it->second == 0)
            s_add *= ~x_vars[it->first];
        else
            s_add *= x_vars[it->first];
    }
    tmp = tmp * s_add;
    tmp.UpdatePathInfo(2, numQubits);
    // tmp.PrintPathInfo();
    return tmp.GetProbability(2, numQubits);
}

std::string BDDQuantumCircuit::Measure() 
{
    ADD tmp = stateVector.SquareTerminalValues();
    tmp.UpdatePathInfo(2, numQubits);
    return tmp.SamplePath(numQubits, 2, "").substr(0, numQubits); 
}

unsigned long long int BDDQuantumCircuit::GetPathCount(long double prob)
{
    ADD tmp = stateVector.SquareTerminalValues();
    tmp.UpdatePathInfo(2, numQubits);
    return tmp.GetPathCount(numQubits, 2, prob); 
}