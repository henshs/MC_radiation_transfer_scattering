// Problem: Monte Carlo simulation of radiative transport via isotropic scattering in a homogeneous medium

#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <string>

using namespace std;

struct quanta{
    double x_coord;
    double y_coord;
    double z_coord;
    double energy;
    double time;
    double radius;
    
    // Constructor to initialize members
    quanta() : x_coord(0.0), y_coord(0.0), z_coord(0.0), energy(1.0), time(0.0), radius(0.0){ }
};

double get_random(double a, double b){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> rand(a, b);
    return rand(gen);
}

int main() {
    // create packets at the center
    // t, x, y, z coordinate
    // propagate
    // track
    // scattering
    //packet array
    //function of index, x, y, z, t, energy
    
    int n_packt = 10000;
    int t_step = 10000;
    
    std::vector<std::vector<quanta>> packet_data(t_step, std::vector<quanta>(n_packt));

    double rand1, rand2, mu, phi, sin_theta;
    double dx_m, dy_m, dz_m;
    
    //move each packet: first step
    for(int i = 0; i < n_packt; i++){
        
        rand1 = get_random(0.0, 1.0);
        rand2 = get_random(0.0, 1.0);
    
        mu = 2 * rand1 - 1; // [-1,-1]
        phi = 2.0 * M_PI * rand2; // [0, 2pi]
        sin_theta = std::sqrt(1.0 - mu*mu);
        
        dx_m = sin_theta * cos(phi);
        dy_m = sin_theta * sin(phi);
        dz_m = mu;
        
        packet_data[0][i].x_coord += dx_m;
        packet_data[0][i].y_coord += dy_m;
        packet_data[0][i].z_coord += dz_m;
        
    }
    
    //move packets
    double rand3, sigma, kmu;
    double rho, tau, mfp, r;
    
    sigma = 1.0; // scattering opacity
    kmu = 0.0; //no absorption
    rho = 1.0; // density
    
    for (int i = 1; i < t_step; i++){
        for(int j = 0; j < n_packt; j++){
            
            rand1 = get_random(0.0, 1.0);
            rand2 = get_random(0.0, 1.0);
            rand3 = get_random(0.0, 1.0);
            
            tau = - log(rand3); // optical depth
            mfp = tau / ((kmu + sigma) * rho); // mean free path
            
            mu = 2 * rand1 - 1; // [-1,-1]
            phi = 2.0 * M_PI * rand2; // [0, 2pi]
            sin_theta = sqrt(1.0 - mu*mu);
            
            dx_m = mfp*sin_theta * cos(phi);
            dy_m = mfp*sin_theta * sin(phi);
            dz_m = mfp*mu;
            
            packet_data[i][j].x_coord = packet_data[i-1][j].x_coord + dx_m;
            packet_data[i][j].y_coord = packet_data[i-1][j].y_coord + dy_m;
            packet_data[i][j].z_coord = packet_data[i-1][j].z_coord + dz_m;
            
            r = sqrt(packet_data[i][j].x_coord*packet_data[i][j].x_coord + packet_data[i][j].y_coord*packet_data[i][j].y_coord + packet_data[i][j].z_coord*packet_data[i][j].z_coord);
            
            packet_data[i][j].radius = r;
    
        }
    }
    
    std::ofstream file("radius_histo.txt");
    
    double rad_max = 0.0;
    double rad_min = 1e6;
    
    for (int j = 0; j < n_packt; j++){
        if (packet_data[t_step - 1][j].radius > rad_max){
            rad_max = packet_data[t_step - 1][j].radius;
        }
        
        if (packet_data[t_step - 1][j].radius < rad_min){
            rad_min = packet_data[t_step - 1][j].radius;
        }
    }
    
    //binning
    int nbins = 100;
    
    std::vector<int> histogram(nbins, 0);
    
    double bin_width = (rad_max - rad_min) / nbins;
    
    for (int j = 0; j < n_packt; j++){
        
        double r = packet_data[t_step - 1][j].radius;
        
        int bin_index = std::min(static_cast<int>((r - rad_min) / bin_width), nbins -1);
        
        histogram[bin_index]++;
    
    }
    
    double rad_mid;
    int histo_debug=0;
    
    for (int i = 0; i < nbins; i++){
        rad_mid = rad_min + (i + 0.5) * bin_width;
        
        file << rad_mid << "   " << histogram[i] << "\n";
        
        histo_debug += histogram[i];
    }
    
    file.close();
    
    cout << "radmax  " << rad_max << endl;
    cout << "radmin  " << rad_min << endl;
    
    cout << "histo debug" << histo_debug << endl;
    
    return 0;
};




//std::ifstream file_rad("radius_histo.txt");
//std::string line;
//std::getline(file_rad, line);
//
//
//while(getline(file_rad, line)){
//    
//}
//
//file_rad.close();





//    for (int i = 0; i < t_step; i++){
//        for(int j = 0; j < n_packt; j++){
//
//            cout << "packet_data[i][j].x_coord" << packet_data[i][j].x_coord << endl;
//            cout << "packet_data[i][j].y_coord" << packet_data[i][j].y_coord << endl;
//            cout << "packet_data[i][j].z_coord" << packet_data[i][j].z_coord << endl;
//            //cout << "packet_data[i][j].time" << packet_data[i][j].time << endl;
//            cout << endl;
//
//        }
//    }


//cout << mfp << endl;

//direction cosine after scattering
//double mu = -1.0 + 2.0 * rand4;

//double rand = get_random(0.0, 1.0);

//cout << rand << endl;

//cout << packet_data[3][3].x_coord << endl;


//    for(int i = 0; i < n_packt; i++){
//        packet_data[i].x_coord = 0.0;
//        packet_data[i].y_coord = 0.0;
//        packet_data[i].z_coord = 0.0;
//        packet_data[i].time = 0.0;
//    };
    
    //cout << "packet_data[1][4].x_coord" << packet_data[1][4].x_coord <<endl;

//    double xmin = -50.0;
//    double xmax = 50.0;
//    double ymin = -50.0;
//    double ymax = 50.0;
//    double zmin = -50.0;
//    double zmax = 50.0;
//
//    int nx = 100;
//    int ny = 100;
//    int nz = 100;
//
//    double dx = (xmax - xmin) / nx;
//    double dy = (ymax - ymin) / ny;
//    double dz = (zmax - zmin) / nz;

    //move packet
//    std::random_device rd;
//    std::mt19937 gen(rd());
//    std::uniform_real_distribution<> random1(0.0, 1.0);
//    std::uniform_real_distribution<> random2(0.0, 1.0);
//    std::uniform_real_distribution<> random3(0.0, 1.0);
//    std::uniform_real_distribution<> random4(0.0, 1.0);

//double t = 0.0;

//cout << dx << dy << dz << endl;


//cout << packet_data[1][0].x_coord << endl;
//cout << sin_theta << endl;
