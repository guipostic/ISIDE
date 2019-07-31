/****************************************************************************
 *   Copyright (C) 2017   Romain Deniau                                     *
 *   Copyright (C) 2017   Jean-Christophe Gelly                             *
 *   Copyright (C) 2017   Yassine Ghouzam                                   *
 *   Copyright (C) 2017   Guillaume Postic                                  *
 *   Copyright (C) 2017   Hubert Santuz                                     *
 *   guillaume.postic@univ-paris-diderot.fr                                 *
 *   jean-christophe.gelly@univ-paris-diderot.fr                            *
 *                                                                          *
 *   This program is free software: you can redistribute it and/or modify   *
 *   it under the terms of the GNU General Public License as published by   *
 *   the Free Software Foundation, either version 3 of the License, or      *
 *   (at your option) any later version.                                    *
 *                                                                          *
 *   This program is distributed in the hope that it will be useful,        *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *   GNU General Public License for more details.                           *
 *                                                                          *
 *   You should have received a copy of the GNU General Public License      *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                          *
 ***************************************************************************/

/* ppot.cpp */

#include "ppot.h"
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

Ppot::Ppot(const string &filename, float _w){
    
    weight=_w;

    // file opening
    ifstream file(filename.c_str(), ios::in);

    if(file){

        // variables
        string aa1,atom1,aa2,atom2;
        vector<float> energies((NUMBER_EVALUES)*NUMBER_INTERPOLATIONS+1);
        string key="",line="";
        istringstream iss;
        float a=0,b=0;
        int j=0;

        // values in Ppot file start at a distance of 0.25 Å (not 0)
        // from 0 to 0.25: energy = 10
        unsigned int start=(NUMBER_INTERPOLATIONS)/2;

        // line by line
        while(getline(file,line)){
            // cast the string into stream (easier to parse)
            iss.str(line);
            // reading of the stream and get back the data
            iss >> aa1 >> atom1 >> aa2 >> atom2;
            // index 0 not used; starts at 1
            //La case 0 aura la meme valeur que la 1ere
            for(unsigned int i=start+1; i<energies.size(); i++)
                if (i%(NUMBER_INTERPOLATIONS)==start || i==(start+1))
                        iss >> energies[i];

            // from 0 to 0.25
            for (unsigned int i=0; i<=start; i++)
                energies[i] = energies[start+1];

            // linear interpolation
            for(unsigned int i=start+1; i<energies.size(); i++){
                if (i%(NUMBER_INTERPOLATIONS)==start  || i==(start+1)){
                    //Derniere case est a 289
                    i==(start+1)?j=start+NUMBER_INTERPOLATIONS:j=i+NUMBER_INTERPOLATIONS;
                    //slope
                    a = (energies[j]-energies[i])/(j-i);
                    //constant term
                    b = energies[i] - i*a;
                }
                else{
                    energies[i] = a*i+b;
                }

            }

            // create object
            Id tmp(aa1,atom1,aa2,atom2,energies);

            // create key
            key=aa1+atom1+aa2+atom2;

            // add to the map
            tab[key]=tmp;

            // flush the stream
            iss.clear();
        }

        // close the file
        file.close();
    }
    else{
        cerr << "Can't open the file named " << filename << endl;
    }
}


const Id &Ppot::operator[](const std::string &key) const{
    return tab.at(key);
}


void Ppot::show(ostream &out) {
    out << "File Ppot: " << endl;
    unordered_map<string, Id>::iterator it = tab.begin();

    for (; it!=tab.end(); ++it)
        out << it->first  << " = "<< it->second << endl;
}


double Ppot::getEnergy(const Id &i, double dist) const{
    double dist2 = dist*PRECISION;
    int access = (int)dist2; // cast is faster than floor()
    double energy=i.getEnergy(access);
    if (energy>2) energy=2; // Limitation du potentiel à 2.
    return (energy*weight);
}


ostream &operator<<( ostream &out, Ppot &d ){
    d.show(out);
    return out;
}
