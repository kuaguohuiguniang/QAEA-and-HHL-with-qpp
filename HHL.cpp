#include "HHLSolver.hpp"
#include "HHLinput.hpp"

// Example usage of HHLSolver
int main() {
    InputHandler input_handler;
    LinearSystem sys = input_handler.load_from_file("data/linear_system.txt");

    int precision = 6;              // Number of clock qubits
    double evolution_time = 2*pi;   // t0
    double rotation_const = 0.25;   // C (must be <= smallest eigenvalue)

    HHLSolver solver(precision, evolution_time, rotation_const);

    qpp::ket solution = solver.solve(sys);
    
    std::cout << "Solution State |x>:\n";
    std::cout << qpp::disp(solution) << std::endl;

    return 0;
}