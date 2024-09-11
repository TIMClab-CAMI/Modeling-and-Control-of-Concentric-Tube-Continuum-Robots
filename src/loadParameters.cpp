#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

#include "loadParameters.h"

using namespace Eigen;

int loadParameters (const std::string &path, std::vector<parameters> &vParameters){
    vParameters.clear();
    std::ifstream indata;
    indata.open(path);
    if(!indata){
        std::cout << "loadParameters()>> Failed to open parameter file located at \"" << path << "\"." << std::endl;
        return -1;
    }
    std::string line;
    int iRow = 0;
    while (std::getline(indata, line)) {
        if(iRow > 0){
            parameters p;
            int iCol = 0;
            std::stringstream lineStream(line);
            std::string cell;
            while (std::getline(lineStream, cell, ',')){
                try{
                    switch(iCol){
                        case 0:
                            p.id = std::stoi(cell);
                            break;
                        case 1:
                            p.name = cell;
                            break;
                        case 2:
                            p.comment = cell;
                            break;
                        case 3:
                            p.E(0) = std::stod(cell);
                            break;
                        case 4:
                            p.E(1) = std::stod(cell);
                            break;
                        case 5:
                            p.E(2) = std::stod(cell);
                            break;
                        case 6:
                            p.mu(0) = std::stod(cell);
                            break;
                        case 7:
                            p.mu(1) = std::stod(cell);
                            break;
                        case 8:
                            p.mu(2) = std::stod(cell);
                            break;
                        case 9:
                            p.Ux(0) = std::stod(cell);
                            break;
                        case 10:
                            p.Ux(1) = std::stod(cell);
                            break;
                        case 11:
                            p.Ux(2) = std::stod(cell);
                            break;
                        case 12:
                            p.l(0) = std::stod(cell);
                            break;
                        case 13:
                            p.l(1) = std::stod(cell);
                            break;
                        case 14:
                            p.l(2) = std::stod(cell);
                            break;
                        case 15:
                            p.l_k(0) = std::stod(cell);
                            break;
                        case 16:
                            p.l_k(1) = std::stod(cell);
                            break;
                        case 17:
                            p.l_k(2) = std::stod(cell);
                            break;
                        case 18:
                            p.rIn(0) = std::stod(cell);
                            break;
                        case 19:
                            p.rIn(1) = std::stod(cell);
                            break;
                        case 20:
                            p.rIn(2) = std::stod(cell);
                            break;
                        case 21:
                            p.rOut(0) = std::stod(cell);
                            break;
                        case 22:
                            p.rOut(1) = std::stod(cell);
                            break;
                        case 23:
                            p.rOut(2) = std::stod(cell);
                            break;
                        case 24:
                            p.Tref(0) = std::stod(cell);
                            break;
                        case 25:
                            p.Tref(1) = std::stod(cell);
                            break;
                        case 26:
                            p.Tref(2) = std::stod(cell);
                            break;
                        case 27:
                            p.Tref(3) = std::stod(cell);
                            break;
                        case 28:
                            p.Tref(4) = std::stod(cell);
                            break;
                        case 29:
                            p.Tref(5) = std::stod(cell);
                            break;
                        case 30:
                            p.Ttip(0) = std::stod(cell);
                            break;
                        case 31:
                            p.Ttip(1) = std::stod(cell);
                            break;
                        case 32:
                            p.Ttip(2) = std::stod(cell);
                            break;
                        case 33:
                            p.offset(0) = std::stod(cell);
                            break;
                        case 34:
                            p.offset(1) = std::stod(cell);
                            break;
                        case 35:
                            p.offset(2) = std::stod(cell);
                            break;
                        case 36:
                            p.offset(3) = std::stod(cell);
                            break;
                        case 37:
                            p.offset(4) = std::stod(cell);
                            break;
                        case 38:
                            p.offset(5) = std::stod(cell);
                            break;
                        case 39:
                            p.w(0) = std::stod(cell);
                            break;
                        case 40:
                            p.w(1) = std::stod(cell);
                            break;
                        case 41:
                            p.w(2) = std::stod(cell);
                            break;
                        case 42:
                            p.w(3) = std::stod(cell);
                            break;
                        case 43:
                            p.w(4) = std::stod(cell);
                            break;
                        case 44:
                            p.w(5) = std::stod(cell);
                            break;
                        case 45:
                            p.w(6) = std::stod(cell);
                            break;
                        case 46:
                            p.w(7) = std::stod(cell);
                            break;
                        case 47:
                            p.lambda = std::stod(cell);
                            break;
                        case 48:
                            p.trajSpeedMax = std::stod(cell);
                            break;
                        case 49:
                            p.trajAccMax = std::stod(cell);
                            break;
                        case 50:
                            p.force(0) = std::stod(cell);
                            break;
                        case 51:
                            p.force(1) = std::stod(cell);
                            break;
                        case 52:
                            p.force(2) = std::stod(cell);
                            break;
                        case 53:
                            p.noiseStd = std::stod(cell);
                            break;
                    }
                }
                catch(...){
                    std::cout << "loadParameters()>> Catch error while parsing ." << std::endl;
                    return -1;
                }
                iCol++;
            }
            if(iCol != 54){
                std::cout << "loadParameters()>> Wrong number of column : " << iCol << " for row " << iRow  << "." << std::endl;
                return -1;
            }
            vParameters.push_back(p);
        }
        ++iRow;
    }
    return 0;
}
void printParameter(parameters &p){
    std::cout << p.id << std::endl;
std::cout << "p.name = " << p.name << std::endl;
std::cout << "p.comment = " << p.comment << std::endl;
std::cout << "p.E1 = " << p.E(0) << std::endl;
std::cout << "p.E2 = " << p.E(1) << std::endl;
std::cout << "p.E3 = " << p.E(2) << std::endl;
std::cout << "p.mu1 = " << p.mu(0) << std::endl;
std::cout << "p.mu2 = " << p.mu(1) << std::endl;
std::cout << "p.mu3 = " << p.mu(2) << std::endl;
std::cout << "p.Ux1 = " << p.Ux(0) << std::endl;
std::cout << "p.Ux2 = " << p.Ux(1) << std::endl;
std::cout << "p.Ux3 = " << p.Ux(2) << std::endl;
std::cout << "p.l1 = " << p.l(0) << std::endl;
std::cout << "p.l2 = " << p.l(1) << std::endl;
std::cout << "p.l3 = " << p.l(2) << std::endl;
std::cout << "p.l_k1 = " << p.l_k(0) << std::endl;
std::cout << "p.l_k2 = " << p.l_k(1) << std::endl;
std::cout << "p.l_k3 = " << p.l_k(2) << std::endl;
std::cout << "p.rIn1 = " << p.rIn(0) << std::endl;
std::cout << "p.rIn2 = " << p.rIn(1) << std::endl;
std::cout << "p.rIn3 = " << p.rIn(2) << std::endl;
std::cout << "p.rOut1 = " << p.rOut(0) << std::endl;
std::cout << "p.rOut2 = " << p.rOut(1) << std::endl;
std::cout << "p.rOut3 = " << p.rOut(2) << std::endl;
std::cout << "p.Trefx = " <<  p.Tref(0) << std::endl;
std::cout << "p.Trefy = " <<  p.Tref(1) << std::endl;
std::cout << "p.Trefz = " <<  p.Tref(2) << std::endl;
std::cout << "p.Trefa1 = " << p.Tref(3) << std::endl;
std::cout << "p.Trefa2 = " << p.Tref(4) << std::endl;
std::cout << "p.Trefa3 = " << p.Tref(5) << std::endl;
std::cout << "p.Ttipx = " << p.Ttip(0) << std::endl;
std::cout << "p.Ttipy = " << p.Ttip(1) << std::endl;
std::cout << "p.Ttipz = " << p.Ttip(2) << std::endl;
std::cout << "p.offsetBeta1 = " <<  p.offset(0) << std::endl;
std::cout << "p.offsetBeta2 = " <<  p.offset(1) << std::endl;
std::cout << "p.offsetBeta3 = " <<  p.offset(2) << std::endl;
std::cout << "p.offsetAlpha1 = " << p.offset(3) << std::endl;
std::cout << "p.offsetAlpha2 = " << p.offset(4) << std::endl;
std::cout << "p.offsetAlpha3 = " << p.offset(5) << std::endl;
std::cout << "p.w0t = " << p.w(0) << std::endl;
std::cout << "p.w0r = " << p.w(1) << std::endl;
std::cout << "p.w1t = " << p.w(2) << std::endl;
std::cout << "p.w1r = " << p.w(3) << std::endl;
std::cout << "p.w2t = " << p.w(4) << std::endl;
std::cout << "p.w2r = " << p.w(5) << std::endl;
std::cout << "p.w3t = " << p.w(6) << std::endl;
std::cout << "p.w3r = " << p.w(7) << std::endl;
std::cout << "p.lambda = " << p.lambda << std::endl;
std::cout << "p.trajSpeedMax = " << p.trajSpeedMax << std::endl;
std::cout << "p.trajAccMax = " << p.trajAccMax << std::endl;
std::cout << "p.forceX = " << p.force(0) << std::endl;
std::cout << "p.forceY = " << p.force(1) << std::endl;
std::cout << "p.forceZ = " << p.force(2) << std::endl;
std::cout << "p.noiseStd = " << p.noiseStd << std::endl;
std::cout << "(width = " << p.rOut-p.rIn << ")" << std::endl;
std::cout << std::endl << std::endl;
}
