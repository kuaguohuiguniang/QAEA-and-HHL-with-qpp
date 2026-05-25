#include "HHLSolver.hpp"
#include "HHLinput.hpp"

/**
 * @brief Runs the example HHL solve using the bundled linear-system input file.
 *
 * The demo loads data/linear_system.txt, configures a small simulated HHL
 * circuit, and prints the postselected solution state.
 */
int main() {
    InputHandler input_handler;
    LinearSystem sys = input_handler.load_from_file("data/linear_system.txt");

    int precision = 6;              ///< Number of clock qubits.
    double evolution_time = 2*pi;   ///< Hamiltonian evolution time t0.
    double rotation_const = 0.25;   ///< Rotation constant C.

    HHLSolver solver(precision, evolution_time, rotation_const);

    qpp::ket solution = solver.solve(sys);
    
    std::cout << "Solution State |x>:\n";
    std::cout << qpp::disp(solution) << std::endl;

    return 0;
}
