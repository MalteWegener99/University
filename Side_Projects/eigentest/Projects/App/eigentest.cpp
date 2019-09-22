#include <Eigen/KroneckerProduct>
#include <Eigen/Sparse>
#include <vector>
#include <iostream>
#include <functional>
#include <math.h>
#include <fstream>

typedef Eigen::SparseMatrix<double> mat;
typedef Eigen::MatrixXd mat2;
typedef std::vector<std::function<double(double)>> function_array;


double source(double x, double y){
    return 20*sin(M_PI * y)*sin(1.5*M_PI*x+M_PI);
}


mat source_flattened(int Nx, int Ny, double h, std::function<double(double, double)> f){
    auto index = [](int nx, int ny, int Nx){return nx+(Nx-1)*ny;};
    int size = (Nx-1)*(Ny-1);
    mat vec(size, 1);
    std::vector<Eigen::Triplet<double>> vals((Nx-1)*(Nx-1));
    int count = 0;
    for(int j = 1; j < Ny; j++){
        for(int i = 1; i < Nx; i++){
            vals[count] = Eigen::Triplet<double>(count, 0, f(i*h, j*h));
            count++;
        }
    }
    vals.shrink_to_fit();
    vec.setFromTriplets(vals.begin(), vals.end());
    return vec;
}

mat make_eigen(int N, int index,double h){
    mat vec(N, 1);
    vec.insert(index, 0) = 1/h/h;
    return vec;
}

mat fill_from_function(int N, double h, std::function<double(double)> f){
    mat vec(N,1);
    for(int i = 1; i <= N; i++){
        vec.insert(i-1,0) = f(i*h);
    }
    return vec;
}

mat make_boundary_vec(int Nx, int Ny, double h, function_array boundaries){
    mat Bx_lower(Nx-1, 1);
    mat Bx_upper(Nx-1, 1);
    mat By_left(Ny-1, 1);
    mat By_right(Ny-1, 1);

    Bx_lower = fill_from_function(Nx-1, h, boundaries[0]);
    By_right = fill_from_function(Ny-1, h, boundaries[1]);
    Bx_upper = fill_from_function(Nx-1, h, boundaries[2]);
    By_left = fill_from_function(Ny-1, h, boundaries[3]);

    return Eigen::kroneckerProduct(make_eigen(Ny-1, 0, h), Bx_lower).eval()+\
           Eigen::kroneckerProduct(make_eigen(Ny-1, Ny-2, h), Bx_upper).eval()+\
           Eigen::kroneckerProduct(By_left, make_eigen(Nx-1, 0, h)).eval()+\
           Eigen::kroneckerProduct(By_right, make_eigen(Nx-1, Nx-2, h)).eval();
}


std::vector<Eigen::Triplet<double>> make_1D_Laplacian(int N, double h){
    std::vector<Eigen::Triplet<double>> vals(N*3);
    int counter = 0;

    for(int i = 0; i < N; i++){
        vals[counter] = Eigen::Triplet<double>(i, i, 2.f/h/h);
        counter++;
        if(i != 0){
            vals[counter] = Eigen::Triplet<double>(i-1, i, -1.f/h/h);
            counter++;
            vals[counter] = Eigen::Triplet<double>(i, i-1, -1.f/h/h);
            counter++;
        }
    }
    vals.shrink_to_fit();
    return vals;
}

mat build_Matrix(int Nx, int Ny, double h){
    mat Dx(Nx-1, Nx-1);
    mat Dy(Ny-1, Ny-1);

    auto Dx_vals = make_1D_Laplacian(Nx-1, h);
    Dx.setFromTriplets(Dx_vals.begin(), Dx_vals.end());

    auto Dy_vals = make_1D_Laplacian(Ny-1, h);
    Dy.setFromTriplets(Dy_vals.begin(), Dy_vals.end());

    mat Iy(Ny-1,Ny-1);
    Iy.setIdentity();

    mat Ix(Nx-1,Nx-1);
    Ix.setIdentity();

    auto L1 = Eigen::kroneckerProduct(Iy, Dx).eval();
    auto L2 = Eigen::kroneckerProduct(Dy, Ix).eval();

    int size = (Nx-1)*(Ny-1);
    return(L1+L2);
}

int normalize(double val, double max, double min){
    double x = 255*((val-min)/(max-min));
    return std::round(x);
}

int main(int argc, char* argv[]){
    double xmax = 2;
    double ymax = 1;
    double h = atof(argv[2]);
    int Nx = (int)(xmax/h);
    int Ny = (int)(ymax/h);
    std::cout << Nx << ", " << Ny << std::endl;
    auto L = build_Matrix(Nx, Ny, h);
    auto source_vec = source_flattened(Nx, Ny, h, std::bind<double>(source, std::placeholders::_1, std::placeholders::_2));
    
    function_array fns(4);
    fns[0] = std::bind<double>([](double x){return sin(0.5*M_PI*x);}, std::placeholders::_1);
    fns[1] = std::bind<double>([](double x){return sin(2*M_PI*x);}, std::placeholders::_1);
    fns[2] = std::bind<double>([](double x){return 0;}, std::placeholders::_1);
    fns[3] = std::bind<double>([](double x){return sin(2*M_PI*x);}, std::placeholders::_1);

    auto boundary_conds = make_boundary_vec(Nx, Ny, h, fns);

    Eigen::SimplicialLDLT<mat> solver;
    L.makeCompressed();
    solver.analyzePattern(L);
    solver.factorize(L);
    std::cout << source_vec.innerSize() << std::endl;
    std::cout << boundary_conds.innerSize() << std::endl;
    std::cout << L.outerSize() << std::endl;
    auto un = solver.solve(source_vec+boundary_conds);

    char* filename;
    if(argc > 1){
        filename = argv[1];
    }
    else{
        filename = new char[20];
        strcpy(filename, "default.pgm");
    }

    auto u = un.toDense();
    double max_val = u.maxCoeff();
    double min_val = u.minCoeff();
    std::cout << max_val << ", " << min_val << std::endl;

    auto norm = std::bind<int>(normalize, std::placeholders::_1, max_val, min_val);

    std::fstream fs;
    fs.open("default.pgm", std::fstream::out | std::fstream::trunc | std::fstream::in);
    std::cout << fs.is_open() << std::endl;
    //print header to file
    fs << "P2\n"<<(Nx-1)<<" "<<(Ny-1)<<"\n255\n";
    for(int j = (Ny-2); j >= 0; j--){
        for(int i = 0; i < (Nx-1); i++){
            fs << norm(u.coeff(i+j*(Nx-1),0)) << " ";
        }
        fs << "\n";

    }
    fs.close();
}