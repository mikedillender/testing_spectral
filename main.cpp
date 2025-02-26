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

const uint16_t width=300;
const double rad_lifetime=1;
double rad_rate = 1/rad_lifetime;
double FRET_denom= 4*(.023*.023);// = 4 sigma^2
double delta_ss=.038;
double FRET_scaling = 30; // (C/dij^6)
double top_density=1;

double num_layers=0;
uint32_t hops=0;


uint16_t energy_resolution=200; // number of pixels on energy axis
double energy_min=2.13;
double energy_span=0.24;
double energy_step=energy_span/energy_resolution;

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
double avg_nr=0;


void save_dots_as_csv(const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
        return;
    }
    int w_edge=5;
    for (uint32_t i=0; i<Np; i++) {
        if(abs(qdx[i])+w_edge>Lbox){i++; continue;}
        if(abs(qdy[i])+w_edge>Lbox){i++; continue;}
        file<<qdz[i]<<",";
        file<<energies[i]<<",";
        file << rates[i][13]<<","<< rates[i][12]<<",";
        for (size_t r = 0; r < 12; ++r) {

            double d_da=(qdx[i]-qdx[j])*(qdx[i]-qdx[j])+(qdy[i]-qdy[j])*(qdy[i]-qdy[j])+(qdz[i]-qdz[j])*(qdz[i]-qdz[j]);
            double rate=rates[i][r]/pow(4/d_da_mag,3);
            file << rate<<","; // Write the value
        }
        file << "\n"; // End the row
    }
    file.close();
    std::cout << "CSV saved to " << filename << std::endl;
}

void save_energy_rates(const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
        return;
    }

    uint32_t E_res=100;
    double en_step=energy_span/E_res;
    vector<vector<double>> Erates=vector<vector<double>>(E_res,vector<double>(E_res,-.0001));
    vector<vector<uint32_t>> points=vector<vector<uint32_t>>(E_res,vector<uint32_t>(E_res,0));
    size_t max_size=12;
    double average_distance=0;
    uint32_t num_distances=0;
    for (uint32_t i=0; i<Np; i++) {
        uint32_t ei=(uint32_t)floor((energies[i]-energy_min)/en_step);
        if(ei>=E_res){continue;}
        size_t num_nns=min(nns[i].size(),max_size);

        for(size_t nn_i=0; nn_i<num_nns; nn_i++){
            //if(j>=nns[id].size()){/*cout<<"not 12 nns\n"; */j=12; continue;}
            uint32_t j=nns[i][nn_i];
            double d_da=(qdx[i]-qdx[j])*(qdx[i]-qdx[j])+(qdy[i]-qdy[j])*(qdy[i]-qdy[j])+(qdz[i]-qdz[j])*(qdz[i]-qdz[j]);
            if(d_da>6){continue;}
            average_distance+=d_da;
            num_distances++;
            uint32_t ej=(uint32_t)floor((energies[j]-energy_min)/en_step);
            if(ej>=E_res){continue;}
            Erates[ei][ej]=(Erates[ei][ej]*points[ei][ej]+rates[i][nn_i])/(points[ei][ej]+1);
            points[ei][ej]++;
        }


    }
    std::cout << "compiled rate matrix, avg dist = "<<average_distance <<", "<<average_distance/num_distances<< std::endl;

    for (uint32_t i=0; i<E_res; i++) {
        for (uint32_t j = 0; j <E_res; j++) {
            file << Erates[i][j]<<","; // Write the value
        }
        file << "\n"; // End the row
    }
    file.close();
    std::cout << "CSV saved to " << filename << std::endl;
}

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
    double energy_mu=2.25; //1240/.55
    double energy_sigma=.030;
    double theta=0;
    double d_z=0;

    vector<double> dix; // unit vector of the dipole * 100
    vector<double> diy; // unit vector of the dipole * 100
    vector<double> diz; // unit vector of the dipole * 100

    random_device rd_e;  // Seed
    mt19937 generator_e(rd_e());  // Random number generator
    normal_distribution<double> e_dist(energy_mu, energy_sigma);
    exponential_distribution<> exp_nr(10/rad_rate);

    // Fill the 3D array with random values from the normal distribution
    for (uint32_t i = 0; i < Np; ++i) {
        energies.push_back(e_dist(generator_e));
        theta = 3.14159*(rand() % 3600)/1800.0;
        d_z = (rand() % 10000 - 5000.0)/5000.0;
        dix.push_back(cos(theta)*sqrt(1-pow(d_z,2)));
        diy.push_back(sin(theta)*sqrt(1-pow(d_z,2)));
        diz.push_back(d_z);
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
        vector<double> rates_i=vector<double>(14,0);

        rates_i[13]=exp_nr(rd_e);
        if(rand()%20==1){
            rates_i[13]=10;
        }
        avg_nr+=rates_i[13];
        rates_i[12]=1+rates_i[13];
        E_d=energies[id];
        r_d[0]=qdx[id]; r_d[1]=qdy[id]; r_d[2]=qdz[id];
        mu_d[0]=dix[id]; mu_d[1]=diy[id]; mu_d[2]=diz[id];
        for(uint32_t j=0;j<12;j++){
            if(j>=nns[id].size()){/*cout<<"not 12 nns\n"; */j=12; continue;}
            uint32_t ia=nns[id][j];
            E_a=energies[ia];
            r_a[0]=qdx[ia]; r_a[1]=qdy[ia]; r_a[2]=qdz[ia];
            mu_a[0]=dix[ia]; mu_a[1]=diy[ia]; mu_a[2]=diz[ia];

            d_da[0]=r_a[0]-r_d[0]; d_da[1]=r_a[1]-r_d[1]; d_da[2]=r_a[2]-r_d[2];
            d_da_mag=d_da[0]*d_da[0]+d_da[1]*d_da[1]+d_da[2]*d_da[2];
            //cout<<"d"<<d_da_mag<<"\n";
            //cout<<"distance";
            kappa= dprod(mu_d,mu_a)-3* dprod(mu_d,d_da)*dprod(mu_a,d_da)/d_da_mag;
            rates_i[j]=FRET_scaling*kappa*kappa*exp(-pow(E_d-E_a-delta_ss,2)/FRET_denom)/pow(.5*(E_d+E_a-delta_ss),4);
            rates_i[j]=rates_i[j]*pow(4/d_da_mag,3);
            rates_i[12]+=rates_i[j];
            //if(d_da_mag<5)
            //    cout<<" - "<<id<<" ("<<E_d<<") -> "<<ia<<" ("<<E_a<<") | rate = "<<rates_i[j]<<", kappa = "<<kappa<<", d = "<<d_da_mag<<"\n";
            //cout<<" | d_da = ("<<d_da[0]<<","<<d_da[1]<<","<<d_da[2]<<") ";
        }
        rates.push_back(rates_i);
        //cout<<id<<", E = "<<E_d<<", total rate = "<<rates[id][12]<<"\n";
    }
    avg_nr=avg_nr/Np;

}

vector<vector<uint32_t>> apd;


uint16_t time_resolution=300; // number of pixels on time axis
double time_max=8;
double time_step=time_max/time_resolution;

random_device rd;
mt19937 rand_gen (rd ());
exponential_distribution<> exp_dist(rad_rate);

void sim_particle(uint32_t i, double t){
    double transfer_time=exp_dist(rd)/rates[i][12]; //should be rad_rate/rates[x][y][z], but rates[][][] is already normalized by base rate
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
    }else if(trans_rand<1+rates[i][13]){
        return;
    }else{
        double cumul_rate=1+rates[i][13];
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
        hops++;
    }
}


void readCSV(const string& filename, vector<double>& x, vector<double>& y, vector<double>& z, bool cut=false, double cut_above=0) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
        return;
    }

    vector<bool> keep;
    vector<uint32_t> new_ind;
    uint32_t kept=0;

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

        if(abs(x_val)>maxx){
            maxx=abs(x_val);
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
        nns.push_back(nearest);
        if(cut){
            new_ind.push_back(kept);
            if(z_val<cut_above){
                x.push_back(x_val);
                y.push_back(y_val);
                z.push_back(z_val);
                kept++;
                keep.push_back(true);
                num_layers=num_layers+z_val;
            } else{
                keep.push_back(false);
            }
        } else{
            x.push_back(x_val);
            y.push_back(y_val);
            z.push_back(z_val);
            num_layers=num_layers+z_val;
        }
        Np++;

    }

    if(cut){
        for(uint32_t i=Np-1; i<Np; i--){
            if(keep[i]){
                for(size_t j=nns[i].size()-1; j<20; j--){
                    if(!keep[nns[i][j]]){
                        nns[i].erase(nns[i].begin()+j);
                    } else{
                        nns[i][j]=new_ind[nns[i][j]];
                    }
                }
            } else{
                nns.erase(nns.begin()+i);
            }
        }
        /*for(uint32_t i=0; i<100; i++) {
            cout<<i<<" : "<<x[i]<<", "<<y[i]<<","<<z[i]<<" : ";
            for(uint32_t j : nns[i]){
                cout<<j<<", ";
            }
            cout<<"\n";
        }*/
        cout<<" kept "<<kept<<"/"<<Np<<" \n";
        Np=kept;
    }
    Lbox=(int32_t)ceil(maxx);
    num_layers=num_layers+Lbox*Np;
    num_layers=(num_layers/Np);
    cout<<"num layers = "<<num_layers<<"\n";
    cout<<"Np = "<<Np<<"\n";
    file.close();
}
int main(/*int argc=0, char** argv=nullptr*/){
    srand((unsigned) time(NULL));
    std::string filename = "dots_z2011_p835_w180_g2_c990.csv"; // Replace with your CSV file path
    //std::string filename = "dots_z1998_p0_w90_g2_c999.csv"; // Replace with your CSV file path
    readCSV(filename, qdx, qdy, qdz, false,-180+3);
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
    int w_edge=10;
    if(width<w_edge*2+2){cout<<"increase width\n"; return 0;}
    cout<<"num layers = "<<num_layers<<"\n";
    uint32_t num_iters=3000000;
    for(uint32_t iter=0; iter<num_iters; iter++){
        uint32_t i = (uint32_t)(rand()%Np);
        if(abs(qdx[i])+w_edge>Lbox){iter--; continue;}
        if(abs(qdy[i])+w_edge>Lbox){iter--; continue;}
        //p=pos((uint16_t)(w_edge+rand()%w_mid),(uint16_t)(w_edge+rand()%w_mid),(uint16_t)(rand()%layers));
        //if(absent[p.x][p.y][p.z]){i--;continue;}
        sim_particle(i,0);
    }
    double avg_hops=((double)hops)/num_iters;
    cout<<"avg hops per excitation : "<<avg_hops<<"\n";
    string name="L"+to_string(Lbox)+"_d"+to_string(int(round(num_layers*100)))+"_C"+to_string(int(FRET_scaling))+"_nr"+to_string(int(round(avg_nr*10)))+".csv";
    string name1="APD_"+name;
    string name2="rates_"+name;
    string name3="E_rates_"+name;
    save_as_csv(apd,name1);
    save_dots_as_csv(name2);
    save_energy_rates(name3);
    return 0;
}