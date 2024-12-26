#include<iostream>
#include<fstream>
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

const uint16_t layers=2;
const uint16_t width=800;
const double lifetime=1;
                          //0 1  2  3 4 5 6 7 8 9 10 11
array<uint16_t,12> rev_dir={9,10,11,8,7,6,5,4,3,0,1, 2};

struct pos{
    pos(uint16_t x,uint16_t y, uint16_t z):x{x},y{y},z{z}{}
    pos(pos &p0):x{p0.x},y{p0.y},z{p0.z}{}
    pos(pos &p0,uint16_t d):x{p0.x},y{p0.y},z{p0.z}{
        go_in_direction(d);
    }
    void go_in_direction(uint16_t d){
        if(d<3){
            if(d==1){ //only change y
                if(z % 2 ==0){ x--;
                }else{ x++; }
            } else if(d==2){ // d==2, change x and y
                if(z % 2 == 0){ // decrease y
                    if (y % 2 == 0){ x--; }
                    y--;
                } else{ // increase y
                    if (y % 2 == 1){ x++; }
                    y++;
                }
            }
            z++;
        }else if(d<9){
            if(d<5){ // d = 3 or 4, up
                if(y%2==0){ x--; }
                if(d==4){ x++; }
                y++;
            }else if(d<7){ // d = 5 or 6, just along x
                if(d==5){ x--;
                }else{ x++; } //if d=6
            }else{ // d = 7 or 8
                if(y%2==0){ x--; }
                if(d==8){ x++; }
                y--;
            }
        }else if(d<12){
            if(d==10){ // only change x
                if(z % 2 ==0){ x--;
                }else{ x++; }
            } else if(d==11){ // d==2, change x and y
                if(z % 2 == 0){ // decrease y
                    if (y % 2 == 0){ x--; }
                    y--;
                } else{ // increase y
                    if (y % 2 == 1){ x++; }
                    y++;
                }
            }
            z--;
        }
    }
    void get_real_pos(double vec[]){
        vec[0]=x;
        vec[1]=y*.86602540378;
        vec[2]=z*0.81649658092;

        if(y%2==1){
            vec[0]+=0.5;
        }
        if(z%2==1){
            vec[0]+=0.5;
            vec[1]-=0.25;
        }
    }
    void print(){cout<<"("<<x<<","<<y<<","<<z<<") ";}
    void print_sub(pos &p1){cout<<"("<<x-p1.x<<","<<y-p1.y<<","<<z-p1.z<<") ";}
    void operator=(const pos& p1){x=p1.x; y=p1.y; z=p1.z;};
    uint16_t x; uint16_t y; uint16_t z;
};


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

double rates[width][width][layers][13];
vector<vector<vector<double>>> energies(width,vector<vector<double>>(width,vector<double>(layers,0)));


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
    signed char dipoles[width][width][layers][3]={}; // unit vector of the dipole * 100

    random_device rd_e;  // Seed
    mt19937 generator_e(rd_e());  // Random number generator
    normal_distribution<double> e_dist(energy_mu, energy_sigma);

    // Fill the 3D array with random values from the normal distribution
    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < width; ++j) {
            for (int k = 0; k < layers; ++k) {
                energies[i][j][k] = e_dist(generator_e);
                theta = 3.14159*(rand() % 360)/180;
                phi = 3.14159*(rand() % 360)/180;
                dipoles[i][j][k][0]=(signed char)(100*cos(theta)*sin(phi));
                dipoles[i][j][k][1]=(signed char)(100*sin(theta)*sin(phi));
                dipoles[i][j][k][2]=(signed char)(100*cos(phi));
            }
        }
    }
    double FRET_denom= 4*(.027*.027);// = 4 sigma^2
    double delta_ss=.038;
    double FRET_scaling = 63; // (C/dij^6)
    double kappa;
    double E_a;
    double E_d;
    double r_d[3]={};
    double r_a[3]={};
    double d_da[3]={};
    double d_da_mag;
    for (uint16_t z=0; z<layers; z++) {
        uint16_t dmin = (z<layers-1)?0:3;
        uint16_t dmax = (z>0)?12:9;
        for (uint16_t x=0; x<width; x++){
            for (uint16_t y=0; y<width; y++) {
                rates[x][y][z][12]=1;
                pos p0=pos(x,y,z);
                p0.get_real_pos(r_d);
                pos p1=pos(p0);
                E_d=energies[x][y][z];
                for (uint16_t d=dmin; d<dmax; d++){
                    p1=pos(p0,d);
                    if(p1.x<width && p1.y<width &&p1.z<layers){
                        p1.get_real_pos(r_a);
                        d_da[0]=r_a[0]-r_d[0]; d_da[1]=r_a[1]-r_d[1]; d_da[2]=r_a[2]-r_d[2];
                        d_da_mag=d_da[0]*d_da[0]+d_da[1]*d_da[1]+d_da[2]*d_da[2];
                        E_a=energies[p1.x][p1.y][p1.z];
                        kappa= dprod(dipoles[x][y][z],dipoles[p1.x][p1.y][p1.z])-3* dprod(dipoles[x][y][z],d_da)*dprod(dipoles[p1.x][p1.y][p1.z],d_da)/d_da_mag;
                        rates[x][y][z][d]=FRET_scaling*kappa*kappa*exp(-pow(E_d-E_a-delta_ss,2)/FRET_denom)/pow(.5*(E_d+E_a-delta_ss),4);
                        rates[x][y][z][12]+=rates[x][y][z][d];
                        //cout<<"\n"; p0.print(); cout<<", "<<E_d<<" -> "; p1.print(); cout<<", "<<E_a<<" | rate = "<<rates[x][y][z][d]<<", kappa = "<<kappa;
                        //cout<<" | d_da = ("<<d_da[0]<<","<<d_da[1]<<","<<d_da[2]<<") ";
                    }
                }
                //cout<<"("<<x<<","<<y<<","<<z<<") E = "<<E_d<<", total rate = "<<rates[x][y][z][12]<<"\n";
            }
        }
    }
}

vector<vector<uint32_t>> apd;

uint16_t energy_resolution=200; // number of pixels on energy axis
double energy_min=2.1;
double energy_span=.3;
double energy_step=energy_span/energy_resolution;

uint16_t time_resolution=200; // number of pixels on time axis
double time_max=5;
double time_step=time_max/time_resolution;

double base_rate = 1/lifetime;
random_device rd;
mt19937 rand_gen (rd ());
exponential_distribution<> exp_dist(base_rate);

void sim_particle(pos p, double t){
    double transfer_time=exp_dist(rd)/rates[p.x][p.y][p.z][12]; //should be base_rate/rates[x][y][z], but rates[][][] is already normalized by base rate
    double trans_rand=((double)rand()/(double)RAND_MAX)*rates[p.x][p.y][p.z][12];// decides how it decays, either FRET or radiative
    //cout<<trans_rand<<"\n";
    if(trans_rand<1){
        uint16_t e_bin=(uint16_t)((energies[p.x][p.y][p.z]-energy_min)/energy_step);
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
            cumul_rate += rates[p.x][p.y][p.z][d];
            if(trans_rand<cumul_rate){
                break;
            }
        }
        /*cout<<"transferred to direction "<<d<<"(rand = "<<trans_rand<<"/"<<rates[p.x][p.y][p.z][12]<<")\n";
        for (uint32_t i=0;i<12;i++){
            cout<<i<<"-"<<rates[p.x][p.y][p.z][i]<<" ";
        }
        cout<<"\n";*/

        p.go_in_direction(d);
        sim_particle(p,t+transfer_time);
    }
}


int main(/*int argc=0, char** argv=nullptr*/){
    setQDs();
    srand((unsigned) time(NULL));


    apd=vector<vector<uint32_t>>(time_resolution,vector<uint32_t>(energy_resolution));
    pos p(0,0,0);
    int w_6=width/6;
    int w_23=width*2/3;
    for(uint32_t i=0; i<1000000; i++){
        p=pos((uint16_t)(w_6+rand()%w_23),(uint16_t)(w_6+rand()%w_23),(uint16_t)(rand()%layers));
        sim_particle(p,0);
    }
    save_as_csv(apd, "output.csv");
    return 0;
}