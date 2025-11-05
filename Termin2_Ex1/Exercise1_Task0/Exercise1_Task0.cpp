#include <iostream>
#include <iomanip>


double f(double x);

double trap(double a, double b, long long N, int num_threads) {
    double h = (b - a) / N;

    double integral = (f(a) + f(b)) / 2;



    #pragma omp parallel for num_threads(num_threads) reduction(+:integral) 
    for (long long i = 1; i < N; i++) {
        double x_i = a + i * h;
        double res = f(x_i);
        integral += res;
    }
    

    integral *= h;

    return integral;
}

int main(int argc, char** argv)
{
    long long samples = 1000000000000000;
    int num_threads = 16;

    if (argc == 3) {
        samples = std::atoi(argv[1]);
        num_threads = std::atoi(argv[2]);
    }

    double result = trap(0, 1, samples, num_threads);

    std::cout << "Approx integral: " << 
        std::setprecision(15) <<
        result << std::endl;
}


double f(double x) {
    return 1.0 / (1.0 + x * x);
}
