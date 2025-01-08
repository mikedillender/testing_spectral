#include<iostream>
#include<fstream>
#include <sstream>
#include<numeric>
#include <random>
#include <vector>
#include <queue>
#include <cstdint>
#include <cmath>
#include <array>
#include <stack>


/*
 * d 0,1,2 are above    |   d 9,10,11 are below
 *                                              y % 2 == 0     y % 2 == 1
 *   ^               3   4                        3  4  -       -  3  4
 *   |     0       5   p   6        9        ->   5  p  6       5  p  6
 *   y   1   2       7   8      10    11          7  8  -       -  7  8
*        x-->
 *           y % 2 = 0                        y % 2 = 1
 *    z % 2 == 0     z % 2 == 1      z % 2 == 0     z % 2 == 1
 *     -  -  -        -  2  -         -  -  -        -  -  2
 *     1  0  -        -  0  1         1  0  -        -  0  1
 *     2  -  -        -  -  -         -  2  -        -  -  -
 */
using namespace std;

const uint16_t layers=3;
const uint16_t width=300;
const double lifetime=1;
double FRET_denom= 4*(.023*.023);// = 4 sigma^2
double delta_ss=.038;
double FRET_scaling = 23; // (C/dij^6)
double top_density=1;


void save_as_csv(const std::vector<std::vector<uint32_t>>& data, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
        return;
    }
    for (const auto& row : data) {
        for (size_t col = 0; col < row.size(); ++col) {
            file << row[col]; // Write the value
            if (col < row.size() - 1) file << ","; // Add a comma except for the last element
        }
        file << "\n"; // End the row
    }
    file.close();
    std::cout << "CSV saved to " << filename << std::endl;
}

vector<double> energies;
vector<vector<double>> rates;

vector<vector<uint32_t>> nns;
uint32_t Np=0;
uint32_t Lbox=0;
std::vector<double> qdx, qdy, qdz;

double dprod(double v1[], double v2[]){
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

double dprod(signed char v1[], double v2[]){
    return (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])/100.0;
}
double dprod(signed char v1[], signed char v2[]){
    return (double)(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])/10000.0;
}

void setQDs(){
    srand((unsigned) time(NULL));
    double energy_mu=2.255; //1240/.55
    double energy_sigma=.030;
    double theta=0;
    double phi=0;

    vector<double> dix; // unit vector of the dipole * 100
    vector<double> diy; // unit vector of the dipole * 100
    vector<double> diz; // unit vector of the dipole * 100

    random_device rd_e;  // Seed
    mt19937 generator_e(rd_e());  // Random number generator
    normal_distribution<double> e_dist(energy_mu, energy_sigma);

    // Fill the 3D array with random values from the normal distribution
    for (uint32_t i = 0; i < Np; ++i) {
        energies.push_back(e_dist(generator_e));
        theta = 3.14159*(rand() % 360)/180;
        phi = 3.14159*(rand() % 360)/180;
        dix.push_back(cos(theta)*sin(phi));
        diy.push_back(sin(theta)*sin(phi));
        diz.push_back(cos(phi));
    }

    double kappa;
    double E_a;
    double E_d;
    double r_d[3]={};
    double r_a[3]={};
    double mu_d[3]={};
    double mu_a[3]={};
    double d_da[3]={};
    double d_da_mag;
    for(uint32_t id=0; id<Np; id++){
        vector<double> rates_i=vector<double>(13,0);
        rates_i[12]=1;
        E_d=energies[id];
        r_d[0]=qdx[id]; r_d[1]=qdy[id]; r_d[2]=qdz[id];
        mu_d[0]=dix[id]; mu_d[1]=dix[id]; mu_d[2]=dix[id];
        for(uint32_t j=0;j<12;j++){
            uint32_t ia=nns[id][j];
            E_a=energies[ia];
            r_a[0]=qdx[ia]; r_a[1]=qdy[ia]; r_a[2]=qdz[ia];
            mu_a[0]=dix[ia]; mu_a[1]=dix[ia]; mu_a[2]=dix[ia];

            d_da[0]=r_a[0]-r_d[0]; d_da[1]=r_a[1]-r_d[1]; d_da[2]=r_a[2]-r_d[2];
            d_da_mag=d_da[0]*d_da[0]+d_da[1]*d_da[1]+d_da[2]*d_da[2];
            //cout<<"d"<<d_da_mag<<"\n";
            //cout<<"distance";
            kappa= dprod(mu_d,mu_a)-3* dprod(mu_d,d_da)*dprod(mu_a,d_da)/d_da_mag;
            rates_i[j]=FRET_scaling*kappa*kappa*exp(-pow(E_d-E_a-delta_ss,2)/FRET_denom)/pow(.5*(E_d+E_a-delta_ss),4);
            rates_i[j]=rates_i[j]*pow(4/d_da_mag,3);
            rates_i[12]+=rates_i[j];
            cout<<id<<" ("<<E_d<<") -> "<<ia<<" ("<<E_a<<") | rate = "<<rates_i[j]<<", kappa = "<<kappa<<", d = "<<d_da_mag<<"\n";
            //cout<<" | d_da = ("<<d_da[0]<<","<<d_da[1]<<","<<d_da[2]<<") ";
        }
        rates.push_back(rates_i);
        //cout<<"("<<x<<","<<y<<","<<z<<") E = "<<E_d<<", total rate = "<<rates[x][y][z][12]<<"\n";
    }

    /*for (uint16_t z=0; z<layers; z++) {
        uint16_t dmin = (z<layers-1)?0:3;
        uint16_t dmax = (z>0)?12:9;
        for (uint16_t x=0; x<width; x++){
            for (uint16_t y=0; y<width; y++) {
                if(absent[x][y][z]){continue;}
                rates[x][y][z][12]=1;
                pos p0=pos(x,y,z);
                p0.get_real_pos(r_d);
                pos p1=pos(p0);
                E_d=energies[x][y][z];
                for (uint16_t d=dmin; d<dmax; d++){
                    p1=pos(p0,d);
                    if(p1.x<width && p1.y<width &&p1.z<layers){
                        if(absent[p1.x][p1.y][p1.z]){continue;}
                        p1.get_real_pos(r_a);
                        d_da[0]=r_a[0]-r_d[0]; d_da[1]=r_a[1]-r_d[1]; d_da[2]=r_a[2]-r_d[2];
                        d_da_mag=d_da[0]*d_da[0]+d_da[1]*d_da[1]+d_da[2]*d_da[2];
                        E_a=energies[p1.x][p1.y][p1.z];
                        kappa= dprod(dipoles[x][y][z],dipoles[p1.x][p1.y][p1.z])-3* dprod(dipoles[x][y][z],d_da)*dprod(dipoles[p1.x][p1.y][p1.z],d_da)/d_da_mag;
                        rates[x][y][z][d]=FRET_scaling*kappa*kappa*exp(-pow(E_d-E_a-delta_ss,2)/FRET_denom)/pow(.5*(E_d+E_a-delta_ss),4);
                        rates[x][y][z][12]+=rates[x][y][z][d];

                    }
                }
            }
        }
    }*/
}

vector<vector<uint32_t>> apd;

uint16_t energy_resolution=200; // number of pixels on energy axis
double energy_min=2.15;
double energy_span=.25;
double energy_step=energy_span/energy_resolution;

uint16_t time_resolution=300; // number of pixels on time axis
double time_max=8;
double time_step=time_max/time_resolution;

double base_rate = 1/lifetime;
random_device rd;
mt19937 rand_gen (rd ());
exponential_distribution<> exp_dist(base_rate);

void sim_particle(uint32_t i, double t){
    double transfer_time=exp_dist(rd)/rates[i][12]; //should be base_rate/rates[x][y][z], but rates[][][] is already normalized by base rate
    double trans_rand=((double)rand()/(double)RAND_MAX)*rates[i][12];// decides how it decays, either FRET or radiative
    //cout<<trans_rand<<"\n";
    if(trans_rand<1){
        uint16_t e_bin=(uint16_t)((energies[i]-energy_min)/energy_step);
        uint16_t t_bin=(uint16_t)((t+transfer_time)/time_step);
        if(e_bin<energy_resolution && t_bin<time_resolution){
            //cout<<"emission at "<<", "<<t_bin<<"\n";
            //cout<<"emission at "<<energies[p.x][p.y][p.z]<<"("<<e_bin<<"), "<<t_bin<<"\n";

            apd[t_bin][e_bin]++;
        }else{
           //cout<<e_bin<<", "<<t_bin<<" out of range\n";
        }
        return;
    }else{
        double cumul_rate=1;
        uint16_t d=0;
        for (; d<12; d++){
            cumul_rate += rates[i][d];
            if(trans_rand<cumul_rate){
                break;
            }
        }
        /*cout<<"transferred to direction "<<d<<"(rand = "<<trans_rand<<"/"<<rates[p.x][p.y][p.z][12]<<")\n";
        for (uint32_t i=0;i<12;i++){
            cout<<i<<"-"<<rates[p.x][p.y][p.z][i]<<" ";
        }
        cout<<"\n";*/


        sim_particle(nns[i][d],t+transfer_time);
    }
}


void readCSV(const string& filename, vector<double>& x, vector<double>& y, vector<double>& z) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
        return;
    }

    std::string line;
    double maxx=0;
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string token;
        vector<uint32_t> nearest;
        double x_val, y_val, z_val;

        // Parse x
        if (std::getline(ss, token, ',')) {
            x_val = std::stod(token);
        } else {
            continue; // Skip if incomplete line
        }

        if(x_val>maxx){
            maxx=x_val;
        }
        // Parse y
        if (std::getline(ss, token, ',')) {
            y_val = std::stod(token);
        } else {
            continue; // Skip if incomplete line
        }

        // Parse z
        if (std::getline(ss, token, ',')) {
            z_val = std::stod(token);
        } else {
            continue; // Skip if incomplete line
        }
        while(std::getline(ss,token,',')){
            nearest.push_back(static_cast<uint32_t>(std::stoul(token)));
        }

        // Add to vectors
        x.push_back(x_val);
        y.push_back(y_val);
        z.push_back(z_val);
        nns.push_back(nearest);
        Np++;
    }
    Lbox=(int32_t)maxx+1;

    file.close();
}

int main(/*int argc=0, char** argv=nullptr*/){
    srand((unsigned) time(NULL));
    std::string filename = "dots_N2880_w30.csv"; // Replace with your CSV file path
    readCSV(filename, qdx, qdy, qdz);
    cout<<" box size is "<<Lbox<<" \n";
    setQDs();

    // Print the vectors for verification
    //std::cout << std::endl;
    /*for (uint32_t i=0;i<Np;i++){
        cout<<i<<" ("<<qdx[i]<<","<<qdy[i]<<","<<qdz[i]<<") : ";
        for(uint32_t j:nns[i]){
            cout<<j<<", ";
        }
        cout<<"\n";
    }*/

    apd=vector<vector<uint32_t>>(time_resolution,vector<uint32_t>(energy_resolution));
    //pos p(0,0,0);
    int w_edge=5;
    if(width<w_edge*2+2){cout<<"increase width\n"; return 0;}
    for(uint32_t iter=0; iter<1000; iter++){
        uint32_t i = (uint32_t)(rand()%Np);
        if(abs(qdx[i])+w_edge>Lbox){iter--; continue;}
        if(abs(qdy[i])+w_edge>Lbox){iter--; continue;}
        //p=pos((uint16_t)(w_edge+rand()%w_mid),(uint16_t)(w_edge+rand()%w_mid),(uint16_t)(rand()%layers));
        //if(absent[p.x][p.y][p.z]){i--;continue;}
        sim_particle(i,0);
    }
    string name="APD_"+to_string(layers)+"L_d"+to_string(int(round(top_density*100)))+"_C"+to_string(int(FRET_scaling))+".csv";
    save_as_csv(apd, name);
    return 0;
}