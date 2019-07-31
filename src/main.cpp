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

#include "PTools/pdbio.h"
#include "PTools/rigidbody.h"
#include "PTools/atomselection.h"
#include "PTools/geometry.h"
#include "PTools/superpose.h"
#include "PTools/basetypes.h"
#include "PTools/rmsd.h"
#include "ppot.h"
#include "solvation.h"
#include "bsheet.h"
#include "alphabet.h"
#include "tabu.h"
#include "distmat.h"
#include "parsetm.h"

#include <algorithm>
#include <iterator>
#include <getopt.h>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <dirent.h>
#include <sys/stat.h>
#include <map>
#include <cctype>

using namespace std;

/* Global const */
const Structural_Alphabet sa("files/blocks");
const Ppot ppot("files/ppot_ncocacb.par", 1.0); // 1.0 = weight
const Solvation sol("files/hydro.par", 32.0); // 32.0 = weight
const Bsheet beta("files/para_x7.mat","files/anti_x7.mat", 2, 9, 100, 0.0); // 3.0 = weight


template<class RandIt>
bool next_k_permutation(RandIt first, RandIt mid, RandIt last){
    typedef typename iterator_traits<RandIt>::value_type value_type;

    sort(mid, last, greater_equal<value_type>());
    return next_permutation(first, last);
}


/* Functor to find the minimum value in a map */
typedef pair<vector<int>, float> MyPairType;
struct CompareSecond{
    bool operator()(const MyPairType &left, const MyPairType &right) const{
        return left.second < right.second;
    }
};
float getMin(map<vector<int>, float> mymap) {
    pair<vector<int>, float> min = *min_element(mymap.begin(), mymap.end(), CompareSecond());
    return min.second;
}
vector<int> getMinKey(map<vector<int>, float> mymap) {
    pair<vector<int>, float> min = *min_element(mymap.begin(), mymap.end(), CompareSecond());
    return min.first; 
}


/* Return the superposition matrix from between the residus start_a-stop_a from the structure 'a' and
 * the residus start_b-stop_b from structure 'b'. */
Matrix get_superposition(const PTools::Rigidbody &a, const PTools::Rigidbody &b, int start_a, int stop_a, int start_b, int stop_b){
    PTools::Rigidbody struct_a = a;
    PTools::Rigidbody struct_b = b;

    /* Select residus for the superposition */
    PTools::AtomSelection ref(struct_a.SelectResRange(start_a, stop_a));
    PTools::AtomSelection mobile(struct_b.SelectResRange(start_b, stop_b));


    /* Discard the oxygen for the superposition */
    ref = ref&(struct_a.SelectAtomType("C")|struct_a.CA()|struct_a.SelectAtomType("N"));
    mobile = mobile&(struct_b.SelectAtomType("C")|struct_b.CA()|struct_b.SelectAtomType("N"));


    /* Create Rigidbody objects */
    PTools::Rigidbody ref_rb(ref.CreateRigid());
    PTools::Rigidbody mobile_rb(mobile.CreateRigid());

    /*Get the superposition matrix */
    return PTools::superpose(ref_rb.Backbone().CreateRigid(), mobile_rb.Backbone().CreateRigid()).matrix;
}


/* Replace a protein block from structure "reference" at the position "position" with the new protein block "new_block"
 * Beware, "position" stands for the central residue of the protein block to switch !!!! */
PTools::Rigidbody switch_PB(PTools::Rigidbody &reference, uint position, PTools::Rigidbody new_block){
    /* Get the structure before the old block
     * "-2" to get the first residue of the block, not the central one */
    PTools::Rigidbody new_struct = reference.SelectResRange(1, position - 2).CreateRigid();

    int last_residue = new_struct.GetAtom(new_struct.Size() - 1).GetResidId();

    /* Superposition of the last residue of "new_struct" and the first residue of the new block */
    Matrix mat = get_superposition(new_struct, new_block, last_residue, last_residue, 1, 1);

    /* Apply the superpose matrix to the new block */
    new_block.ApplyMatrix(mat);
    new_block.Renumber(0);

    /* Merge only if it is not the first block of the sequence to mutate */
    PTools::Rigidbody c;
    if (position - 2 > 1) {
        c = Merge(new_struct, new_block, 1, position - 2, 1, 4);
    }
    else {
        new_block.Renumber(2);
        new_block.Renumber(1);
        c = new_block;
    }
    c.Renumber();

    /* End of the first part of the switch
    Now, re-attach the rest of the reference structure with the block switched*/

    /* Get the structure after the old block
     * "+2" to get the last residue of the block, not the central one */
    PTools::Rigidbody end = reference.SelectResRange(position + 2, 
                                                     reference.GetAtom(reference.Size() - 1).GetResidId())
                                     .CreateRigid();
    end.Renumber();

    last_residue = c.GetAtom(c.Size() - 1).GetResidId();
    mat = get_superposition(c, end, last_residue, last_residue, 1, 1);

    /* Apply the superpose matrix to the rest of the structure */
    end.ApplyMatrix(mat);

    /* Merge the 2 blocks */
    PTools::Rigidbody c2 = Merge(c, end, 1, last_residue - 1, 1, end.Size());
    c2.Renumber();

    return c2;
}


/* Attach 1 block (block2) to a 3D structure of blocks . the superposition before the merge is dictated by the variable "position" of the structure :
 *  - central position (i) or i-1,i-2,i+1,i+2
 *  pos_structure is between 2 and 5 */
PTools::Rigidbody attach_PB(const PTools::Rigidbody &structure, const PTools::Rigidbody &block2, uint pos_structure) {
    PTools::Rigidbody res;
    PTools::Rigidbody tmp = block2;

    /* Position (in AA) for the superposition
      pos_structure is between {2,5} */
    uint pos_block2 = pos_structure - 1;
    uint nbr_PB_remove_structure = Structural_Alphabet::nbr_aa_per_block - pos_structure;

    /* Create Selection to superpose */
    Matrix mat = get_superposition(structure, tmp, structure.Size_residus() - nbr_PB_remove_structure,
                                   structure.Size_residus() - nbr_PB_remove_structure, pos_block2, pos_block2);

    /* Apply the superpose matrix to the 2nd block */
    tmp.ApplyMatrix(mat);

    /* Merge the 2 blocks
     * Take 1 to (nbr_residus- 3) residus from block 1 and the 3 last from block 2 */
    res = Merge(structure, tmp, 1, structure.Size_residus() - nbr_PB_remove_structure - 1,
                pos_block2, Structural_Alphabet::nbr_aa_per_block);

    if (pos_structure == 2)
        res.Renumber_res();
    else
        res.Renumber();

    res.set_block_residue(structure.get_block(structure.Size_residus() - nbr_PB_remove_structure),
                          structure.Size_residus() - nbr_PB_remove_structure);
    res.set_inpos_block_per_residue(pos_structure, structure.Size_residus() - nbr_PB_remove_structure);

    return res;
}


/* Create a 3D structure from the PB sequence (seq) and the structural alphabet (sa)
 * Return the structure */
PTools::Rigidbody Create_structure(const string &seq){
    if (seq.size() == 0 ) {
        return PTools::Rigidbody();
    }

    /* First Block */
    PTools::Rigidbody res = sa.get_block(seq[0]);

    /* For each others blocks */
    for(uint i = 1; i < seq.size(); i++) {
        /* Next block */
        PTools::Rigidbody next = sa.get_block(seq[i]);
        /* Attach it to the structure */
        res = attach_PB(res, next, 3);
    }

    return res;
}


void cart_product(
    vector<vector<char> > &rvvi,  // final result
    vector<char> &rvi,   // current result 
    vector<vector<char> >::const_iterator me, // current input
    vector<vector<char> >::const_iterator end) // final input
{
    if(me == end) {
        // terminal condition of the recursion. We no longer have
        // any input vectors to manipulate. Add the current result (rvi)
        // to the total set of results (rvvvi).
        rvvi.push_back(rvi);
        return;
    }

    // need an easy name for my vector-of-ints
    const vector<char>& mevi = *me;
    for(vector<char>::const_iterator it = mevi.begin(); it != mevi.end(); it++) {
        // final rvi will look like "a, b, c, ME, d, e, f"
        // At the moment, rvi already has "a, b, c"
        rvi.push_back(*it);  // add ME
        cart_product(rvvi, rvi, me+1, end); // add "d, e, f"
        rvi.pop_back(); // clean ME off for next round
    }
}


/* Return a map<vector<int>, float> containing the *available* possibilities for the mutations.
 * The key of the map is a vector of indexes of the tabu list. The length of this vector is dictated by the number of mutations
 * The value of the map is the energy according to the mutations done. */
typedef map<vector<int>, float> Map_index_ener;

Map_index_ener generate_indexes(int size_seq_PB, const Tabu &tabulist, uint nbr_possibilities,
                                const string & locked_pos, bool montecarlo, float threshold, float threshold2){

    /* Now, generate all possibilities between all positions and all blocks */
    Map_index_ener index_energies;

    /* Create a vector with a size equal to the size of the BP sequence */
    vector<int> vseq;
    for(int i = 1; i <= size_seq_PB; ++i)
       if(locked_pos[i-1]=='O')
          vseq.push_back(i);

    /* Generate all permutations from the vseq and the number of mutations
     * [1,2],[1,3],[1,4],...,[10,11],[10,12],... with double mutations
     * Store them in vector of vector of int
     */
    vector<vector<int> > permut_pos;
    do{
        vector<int> tmp(vseq.begin(), vseq.begin() + 2); // XXX 2 for double mutations
        permut_pos.push_back(tmp);
    }
    while(next_k_permutation(vseq.begin(), vseq.begin() + 2, vseq.end())); // XXX 2 for double mutations

    /* Monte Carlo */
    uint limitA = 0;
    uint limitB = permut_pos.size();

    if (montecarlo){
        limitA = rand() % limitB;
        limitB = limitA+1;
    }

    /* For every positions */
    for (uint i = limitA; i < limitB; i++) {
        vector<vector<char> > best_blocks;

        for (uint j = 0; j < permut_pos[i].size(); j++) {
            /* Find the 'nbr_possibilites' best available block with highest proba at position 'permut_pos[i][j]' in tabu list */
            vector<char> tmp = tabulist.find_best_blocks(permut_pos[i][j], nbr_possibilities, threshold, threshold2);
            best_blocks.push_back(tmp);
        }

        /* Generate combinations for theses best blocks */
        vector<vector<char> > out_blocks;
        vector<char> combi_tmp;
        cart_product(out_blocks, combi_tmp, best_blocks.begin(), best_blocks.end());

        /* Now, generate all possibilities between all positions and best blocks */
        for (uint j = 0; j < out_blocks.size(); j++) {

            vector<int> list_index; /*Key of the map */
            for(uint k = 0; k < permut_pos[i].size(); k++) {
                /* Calculate the good index in the tabu list */
                int index_tabu_list = (permut_pos[i][k] - 1)*sa.size() + correspond(out_blocks[j][k]);
                if (index_tabu_list == -1) {                   
                    cout << out_blocks[j][0] << " - " << out_blocks[j][1] << endl;
                    cout << permut_pos[i][k] -1 << " --- " << out_blocks[j][0] << ", " << permut_pos[i].size()  << endl;
                }
                list_index.push_back(index_tabu_list);
            }
            index_energies[list_index] = 0;
        }
    }

    return index_energies;
}


/* Since it is not possible to use Map::iterator with openMP, we need
 * to create a vector of iterator. This vector will be paralelize with openmp
 * Returns a vector of iterator. Take the map object in argument */
vector<Map_index_ener::iterator> create_it_array(Map_index_ener & map_index){
    int it_array_size = map_index.size();
    vector<Map_index_ener::iterator> vec_it(it_array_size);

    Map_index_ener::iterator allo_it = map_index.begin();
    for (int j = 0; j < it_array_size; ++j)
        vec_it[j] = allo_it++;

    return vec_it;
}


void print_information(int i, string new_positions, string new_blocks, const PTools::Rigidbody &new_struct,
                       float new_energy, float new_rmsd, ofstream & printout, ofstream & weblogo,
                       float new_tmscore, float new_gdt){

    printout << i << " "
             << new_positions << " "
             << new_blocks << " "
             << new_struct.get_all_blocks() << " "
             << new_energy << " "
             << new_rmsd << " "
             << new_gdt << " "
             << new_tmscore << endl;

    weblogo << ">" << i << endl << new_struct.get_all_blocks() << endl;

    string separator = "   ";
    string special = (new_energy > -100) ? "    " :
                     (new_energy< -1000) ? "  " : "   ";

    cout << setfill('0') << setw(5) << i << separator
         << setfill(' ') << setw(6) << new_positions << separator
         << setfill(' ') << setw(5) << new_blocks << separator
         << setfill(' ') << setw(new_struct.get_all_blocks().size()) <<  new_struct.get_all_blocks() << special
         << setfill(' ') << setw(6) << fixed << setprecision(2) << new_energy << separator
         << setfill(' ') << setw(8) << fixed << setprecision(2) << new_tmscore << separator
         << setfill(' ') << setw(8) << fixed << setprecision(2) << new_gdt << separator
         << setfill(' ') << setw(6) << fixed << setprecision(2) << new_rmsd << endl;
}


Map_index_ener do_step0(const Tabu &tabulist, PTools::Rigidbody &ref, const vector<string> residus_seq,
    int size_seq_PB, int nbr_best_choices, const string & locked_pos, bool montecarlo, float threshold, float threshold2){

    /* Generate all possibilities and store them in a map */
     Map_index_ener index_energies = generate_indexes(size_seq_PB, tabulist, nbr_best_choices, locked_pos, montecarlo, threshold, threshold2);

    /* Create a vector of map::iterator for parallelization
     * It takes less than 1s */
    vector<Map_index_ener::iterator> it_array = create_it_array(index_energies);

    uint i;

    #pragma omp parallel
    {
        #pragma omp for schedule(static) private(i)
        for (i = 0; i < it_array.size(); i++) { /*For every possibility */
        
            /* XXX filter for max distance between 2 mutations (only between fisrt and second mutation) */
            uint positions[2];
            uint indice=0;
            vector<int>::const_iterator j;
            for(j = (it_array[i]->first).begin(); j != (it_array[i]->first).end(); ++j) {
               positions[indice]=tabulist[*j].get_pos() + 2;
               indice++;
            }

            if(positions[1]-positions[0]<20){ // MAX_MUT_DIST : maximum distance between 2 mutations
                /* Generation structure */
                PTools::Rigidbody tmp_struct = ref;
                vector<int>::const_iterator j;
                for(j = (it_array[i]->first).begin(); j != (it_array[i]->first).end(); ++j) {
                    PTools::Rigidbody block_3D = sa.get_block(tabulist[*j].get_block());
                    tmp_struct = switch_PB(tmp_struct, tabulist[*j].get_pos() + 2, block_3D);
                }
                tmp_struct.change_all_resnames(residus_seq);

                /* Calculate energy */
                it_array[i]->second = tmp_struct.compute_energies(ppot,sol,beta);
            }
        }
    }

    return index_energies;
}


map<char,string> gen_1to3() {
     map<char,string> m;
     m['A'] = "ALA"; m['R'] = "ARG"; m['N'] = "ASN"; m['D'] = "ASP";
     m['C'] = "CYS"; m['E'] = "GLU"; m['Q'] = "GLN"; m['G'] = "GLY";
     m['H'] = "HIS"; m['I'] = "ILE"; m['L'] = "LEU"; m['K'] = "LYS";
     m['M'] = "MET"; m['F'] = "PHE"; m['P'] = "PRO"; m['S'] = "SER";
     m['T'] = "THR"; m['W'] = "TRP"; m['Y'] = "TYR"; m['V'] = "VAL";
     return m;
}


map<string,char> gen_3to1() {
     map<string,char> m;
     m["ALA"] = 'A'; m["ARG"] = 'R'; m["ASN"] = 'N'; m["ASP"] = 'D';
     m["CYS"] = 'C'; m["GLU"] = 'E'; m["GLN"] = 'Q'; m["GLY"] = 'G';
     m["HIS"] = 'H'; m["ILE"] = 'I'; m["LEU"] = 'L'; m["LYS"] = 'K';
     m["MET"] = 'M'; m["PHE"] = 'F'; m["PRO"] = 'P'; m["SER"] = 'S';
     m["THR"] = 'T'; m["TRP"] = 'W'; m["TYR"] = 'Y'; m["VAL"] = 'V';
     return m;
}


vector<string> mkdir_results(string root_rep, bool fix_pos) {
    
    vector<string> dirs;

    if(fix_pos) {
        dirs.push_back(root_rep+"_fixpos"+"_all");
        dirs.push_back(root_rep+"_fixpos"+"_best");
    } else {
        dirs.push_back(root_rep+"_all");
        dirs.push_back(root_rep+"_best");
    }

    for(unsigned int i=0; i< dirs.size(); i++) {
        mkdir(dirs[i].c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }
    
    return dirs;
}


int errwarn(string jobname, string ropt, string aaseq, string pdbname, string nopt, string topt, string file_proba,
             bool hiprobaPB, string yopt, float threshold, float threshold2, uint& nbr_possibilities, bool montecarlo,
             int& size_tabu, string kopt, string dopt, string aopt, string seq_locked_pbe){
    if (jobname.empty()){
        cerr << "* Error: please enter a job name (-j option)" << endl;
        return 1;
    }

    if (!ropt.empty() && ropt.substr(0,2) != "0."){
        cerr << "* Error: -r argument must belong to [0;1[" << endl;
        return 1;
    }

    if (!yopt.empty() && yopt.substr(0,2) != "0."){
        cerr << "* Error: -y argument must belong to [0;1[" << endl;
        return 1;
    }

    if (!aaseq.empty() && !pdbname.empty()){
        cerr << "* Error: -i and -c options are not compatible" << endl;
        return 1;
    }

    if (!seq_locked_pbe.empty() && !yopt.empty()){
        cerr << "* Error: -l and -y options are not compatible" << endl;
        return 1;
    }

    if ( (isdigit(nopt[0]) && (atoi(nopt.c_str()) <= 0 || atoi(nopt.c_str()) > 15)) || !isdigit(nopt[0]) ){
        cerr << "* Error: -n argument must belong to [1;16[" << endl;
        return 1;
    }

    if ( (isdigit(topt[0]) && atoi(topt.c_str()) < 0) || !isdigit(topt[0])){
        cerr << "* Error: -t argument must be a positive integer" << endl;
        return 1;
    }

    if (file_proba.empty() && (hiprobaPB || !ropt.empty())){
        cerr << "* Error: -x and -r options require -f option" << endl;
        return 1;
    }

    if (threshold > 0 && threshold2 > 0){
        cerr << "* Error: -r and -q options are not compatible" << endl;
        return 1;
    }

    if (file_proba.empty()){
        cerr << "* Warning: -f option not used: all the 16 PBs will be tested (will slow down ISIDE)\n" << endl;
        nbr_possibilities = 15;
        return 0;
    }
    else if (nbr_possibilities == 0){
        cerr << "* Error: -n option must be > 0" << endl;
        return 1;
    }

    if (montecarlo && size_tabu > 0){
        cerr << "* Warning: Monte Carlo activated (-m option): tabu list size (-t option) is set to 0\n" << endl;
        size_tabu = 0;
        return 0;
    }

    if (!montecarlo && (!kopt.empty() || !dopt.empty() || !aopt.empty())){
        cerr << "* Warning: -k, -d, and -a arguments (Monte Carlo parameters) are useless without -m option\n" << endl;
        return 0;
    }

    if ((threshold > 0 || threshold2 > 0) && nbr_possibilities < 15){
        cerr << "* Warning: with -r or -q options, \"-n 15\" is recommended\n" << endl;
        return 0;
    }

    return 0;
}


/* (not imortant) parser aimed at making the -l option more intuitive */
int parse_lock(string& seq_locked_pbe, bool die){
    string temp_locked_pbe = seq_locked_pbe;
    uint lastx = 0;
    for (size_t i=3; i<seq_locked_pbe.size(); ++i){ // '<=' useless
        if (seq_locked_pbe[i] == 'X'){
            uint diffx = i-lastx;
            if ( lastx != 0 && (diffx == 2 || diffx == 3) ){ // do nothing or die!
                if (die){
                    cerr << "* Error: 'X's must be separated by 3 or more 'O's" << endl;
                    return 1;
                }
            }
            else{
                temp_locked_pbe[i-3] = 'X';
                temp_locked_pbe[i-2] = 'X';
                temp_locked_pbe[i-1] = 'X';
                temp_locked_pbe[i] = 'O';
                lastx = i;
            }
        }
        else if (seq_locked_pbe[i] == 'O'){
            temp_locked_pbe[i] = 'O';
        }
        else{
            cerr << "* Error: -l arg string must contain only X or O" << endl;
            return 1;
        }
    }
    seq_locked_pbe = temp_locked_pbe;
    return 0;
}


int best_start(string& file_proba, string& seq_assign_pbe, vector<float>& seq_assign_proba, size_t seqlength){

    if (file_proba.empty()){
        cerr << "* Error: -x option requires -f option" << endl;
        return 1;
    }

    if (!seq_assign_pbe.empty()){
        cerr << "* Error: -x and -p options are not compatible" << endl;
        return 1;
    }

    string alphabet[] = {"a", "b", "c", "d", "e","f","g","h","i","j","k","l","m","n","o","p"};
    
    float hival = 0.0;
    float probaPB;
    ifstream inFile;
    
    int i = 0;
    int j = 0;
    uint linecount = 0;
    
    inFile.open(file_proba.c_str());
    while ( inFile >> probaPB ){
        if (probaPB > hival){
            hival = probaPB;
            j = i;
        }

        if (i == 15){
            seq_assign_pbe+=alphabet[j];
            seq_assign_proba.push_back(hival);

            i = 0;
            hival = 0;
            linecount++;
        }
        else{
            i++;
        }
    }
    inFile.close();

    if (linecount != seqlength){
        cerr << "* Error: PB matrix (-f option) and protein sequence must have the same length" << endl;
        cerr << "matrix=" << linecount << "; seq=" << seqlength << endl;
        return 1;
    }
 
    size_t len = seq_assign_pbe.length()-4;
    seq_assign_pbe = seq_assign_pbe.substr(2,len);

    return 0;
}


void lockrules(string& seq_assign_pbe, vector<float>& seq_assign_proba, string& seq_locked_pbe, float& locksse){
    for (size_t i = 0; i < seq_assign_pbe.size(); ++i){
        if (seq_assign_pbe[i] == 'm' || seq_assign_pbe[i] == 'd'){
            if (seq_assign_proba[i] > locksse){
                seq_locked_pbe[i+2] = 'X';
            }
        }
        if (seq_assign_pbe[i] == 'c' || seq_assign_pbe[i] == 'e'){
            if ( i>0 && (i+1)<seq_assign_pbe.size() ){
                if (seq_assign_pbe[i-1] == 'd' || seq_assign_pbe[i+1] == 'd'){
                    //if (seq_assign_proba[i] > locksse){
                    if (seq_assign_proba[i] > locksse &&
                        seq_assign_proba[i-1] > locksse &&
                        seq_assign_proba[i+1] > locksse){

                        seq_locked_pbe[i+2] = 'X';
                    }
                }
            }
        }
        if (seq_assign_pbe[i] == 'l' || seq_assign_pbe[i] == 'n'){
            if ( i>0 && (i+1)<seq_assign_pbe.size() ){
                if (seq_assign_pbe[i-1] == 'm' || seq_assign_pbe[i+1] == 'm'){
                    //if (seq_assign_proba[i] > locksse){
                    if (seq_assign_proba[i] > locksse &&
                        seq_assign_proba[i-1] > locksse &&
                        seq_assign_proba[i+1] > locksse){

                        seq_locked_pbe[i+2] = 'X';
                    }
                }
            }
        }
    }
}


int main(int argc, char *argv[]) {
    srand (time(NULL));

    /* Check if TMscore executable exists */
    ifstream infile("TMscore/TMscore");
    if ( !infile.good() ){
        cerr << "\n* Error: TMscore executable not found" << endl;
        return 1;
    }

    string optlist =
        "                              , ____                               \n"
        "                              ~(__  \\.                             \n"
        "     .--------__               :o \\. |                 __--------. \n"
        "    ( _________ \\              `.  |_|                / _________ )\n"
        "     ( _________ `o-___________-`-.-'-_____________-o' _________ ) \n"
        "      ( _________  `-----------(.v(.)-------------'  _________ )   \n"
        "        \\____________________  `-   -'  ______________________/    \n"
        "           \\___________________'\\ * /`_____________________/       \n"
        "                          ______: . :                              \n"
        "                         / _____ /   )                             \n"
        "                        : /    /'  /'                              \n"
        "                        : l   /  /                                 \n"
        "                        :)   (  :---._                             \n"
        "                     ..//     `--------O\\..                        \n"
        "                                                                 \n"
        "    ISIDE: alphabet-based metaheuristic for structure prediction \n"
        "\n\n"
        "    Usage: ./ISIDE [options]\n\n"
        "    Options:\n"
        "    -i <string>    Input sequence (amino acids)\n"
        "    -t <int>       Tabu list size (default=0)\n"
        "    -n <int>       Number of most probable PBs tested for each position (default=4)\n"
        "    -s <int>       Steps, i.e. the number of iterations (default=500)\n"
        "    -f <file>      File containing the PBs probabilities\n"
        "    -c <file>      Comparison PDB file\n"
        "    -o <file>      Original structure file: start from an existing structure\n"
        "    -p <string>    PBs initial sequence\n"
        "    -l <string>    Locking for every position: sequence of 'X's (locked) and 'O's (open)\n"
        "    -j <string>    Job name\n"
        "    -x             Start from the most probable PB sequence (instead of a sequence of 'd')\n"
        "    -y <float>     Same as -x option, with most probable (>threshold) secondary structures locked\n"
        "    -r <float>     Threshold probability (default=0.0): test only PBs with a higher probability\n"
        "    -q <float>     Cumulated threshold probability (default=0.0); alternative to -r option\n"
        "    -m             Monte Carlo activated (FOR TESTING PURPOSE ONLY!)\n"
        "    -k <float>     Monte Carlo: initial Ck value (default=2.0)\n"
        "    -a <float>     Monte Carlo: Ck decrease factor, i.e. 'cooling' (default=0.99)\n"
        "    -d <int>       Monte Carlo: Ck decreases by [-a] every [-d] steps (default=2e+5)\n"
        "    -h             Help\n"
        "\n";

    string aaseq; // -i opt    

    string file_proba; // -f opt
    stringstream out;
    bool montecarlo = false; // -m opt
    bool hiprobaPB = false; // -x opt

        // -c opt     -o opt        -p opt         -l opt        -j opt
    string pdbname, pb_pdbname, seq_assign_pbe, seq_locked_pbe, jobname;
    string topt, nopt, sopt, kopt, dopt, aopt, ropt, qopt, yopt;

    int opt;
    while ((opt = getopt(argc,argv,"hxi:k:d:a:t:n:s:r:q:f:c:o:p:l:j:h:y:m")) != EOF){
        switch(opt){
            case 'i': aaseq = optarg; break;
            case 'm': montecarlo = true; break;
            case 'x': hiprobaPB = true; break;
            case 'y': yopt = optarg; break;
            case 't': topt = optarg; break;
            case 'n': nopt = optarg; break;
            case 's': sopt = optarg; break;
            case 'k': kopt = optarg; break;
            case 'd': dopt = optarg; break;
            case 'a': aopt = optarg; break;
            case 'r': ropt = optarg; break;
            case 'q': qopt = optarg; break;
            case 'f': file_proba = optarg; break;
            case 'c': pdbname = optarg; break;
            case 'o': pb_pdbname = optarg; break;
            case 'p': seq_assign_pbe = optarg; break;
            case 'l': seq_locked_pbe = optarg; break;
            case 'j': jobname = optarg; break;
            case 'h': fprintf(stderr, "%s", optlist.c_str()); return 0;
        }
    }

    if (argc == 1){
        fprintf(stderr, "%s", optlist.c_str());
        return 0;
    }
    cout << endl;

    uint nbr_possibilities = (!nopt.empty()) ? atoi(nopt.c_str()) : 4;
    int size_tabu = (!topt.empty()) ? atoi(topt.c_str()) : 0;
    int steps = (!sopt.empty()) ? atoi(sopt.c_str()) : 500;
    float threshold = (!ropt.empty()) ? atof(ropt.c_str()) : 0.0;
    float threshold2 = (!qopt.empty()) ? atof(qopt.c_str()) : 0.0;
    float locksse = (!yopt.empty()) ? atof(yopt.c_str()) : 0.0;
    if (locksse > 0){
        hiprobaPB = true;
        cerr << "* Note: -y option activates -x option\n" << endl;
    }

    /* Arguments for the Monte Carlo */
    float ck = (!kopt.empty()) ? atof(kopt.c_str()) : 2.0;
    int decrease =  (!dopt.empty()) ? atoi(dopt.c_str()) : 200000;
    float alpha = (!aopt.empty()) ? atof(aopt.c_str()) : 0.99;

    /* Errors and warnings */
    if ( errwarn(jobname, ropt, aaseq, pdbname, nopt, topt, file_proba, hiprobaPB, yopt, threshold,
                  threshold2, nbr_possibilities, montecarlo, size_tabu, kopt, dopt, aopt, seq_locked_pbe) ){
        return 1;
    }

    string printout_name = string(jobname) + ".txt";
    string weblogo_name = string(jobname) + "_weblogo.fasta";
    ofstream printout(printout_name.c_str());
    ofstream weblogo(weblogo_name.c_str());

    /* 3-letter residue list */
    vector<string> vec_seq_residues;
    if (!aaseq.empty()){
        map<char,string> one_to_three = gen_1to3();

        vector<char> vec_aaseq1(aaseq.begin(), aaseq.end());

        for(size_t i = 0; i<vec_aaseq1.size(); ++i){
            char keyi = vec_aaseq1[i];
            vec_seq_residues.push_back(one_to_three[keyi]);
        }
    }
    else{
        PTools::Rigidbody original(pdbname);
        original.Renumber();
        vec_seq_residues = original.get_all_resnames();
    }


    /* if start from the most probable PBs sequence */
    vector<float> seq_assign_proba; // related to -y option
    string raw_seq_locked_pbe = "n/a";
    if (hiprobaPB){
        if ( best_start(file_proba, seq_assign_pbe, seq_assign_proba, vec_seq_residues.size()) == 1 ){
            map<string,char> three_to_one = gen_3to1();
            string oneletterseq;
            for(size_t i = 0; i<vec_seq_residues.size(); ++i){
                string keyi = vec_seq_residues[i];
                oneletterseq+=three_to_one[keyi];
            }
            cerr << "Seq: " << oneletterseq << endl;

            return 1;
        }
    }

    if (seq_assign_pbe.empty())
        for (size_t i = 1; i <= vec_seq_residues.size()-4; ++i) // -4 = the two flanking "ZZ"
            seq_assign_pbe+="d";

    if (!seq_locked_pbe.empty()){ // if -l option
        raw_seq_locked_pbe = seq_locked_pbe;
        if ( parse_lock(seq_locked_pbe, 1) == 1 ) return 1;
    }
    else{
        //for (size_t i = 1; i <= vec_seq_residues.size()-2; ++i){ // useless
        for (size_t i = 1; i <= vec_seq_residues.size(); ++i){
            seq_locked_pbe+="O";
        }
        if (locksse > 0){ // if -y option
        /* lock regular secondary structure elements */
            lockrules(seq_assign_pbe, seq_assign_proba, seq_locked_pbe, locksse);

            /* In principle, this one should never return 1 */
            raw_seq_locked_pbe = seq_locked_pbe;
            if ( parse_lock(seq_locked_pbe, 0) == 1 ) return 1;
        }
    }


    int size_seq_PB = seq_assign_pbe.size();
    Tabu tabulist = (file_proba.empty()) ? Tabu(size_tabu, size_seq_PB) : Tabu(file_proba, size_tabu);


    /* Build initial structure */
    PTools::Rigidbody Sref = (!pb_pdbname.empty()) ? PTools::Rigidbody(pb_pdbname) : Create_structure(seq_assign_pbe);

    Sref.change_all_resnames(vec_seq_residues);
    Sref.set_all_blocks(seq_assign_pbe);

    /* Create the two results dir */
    string job_dir = "ISIDE_" + jobname;
    vector<string> result_dirs = mkdir_results(job_dir, false);

    Sref.init_blocks_edges();
    cout << "Initial PBs sequence: " << Sref.get_all_blocks() << endl;
    cout << "Locked positions:     " << raw_seq_locked_pbe << endl;
    cout << "\nStart calculations...\n" << endl;

    /* Initialize the matrix */
    tabulist.init(seq_assign_pbe);


    /* Step 0
       Create all the structures and calculate the energy without modifying the initial sequence
       For every position:
       - Create a structure for all the other PBs (or only those with the highest probabilities; see options)
       - Calculate the energies of these structures and save them in the tabu list
       The mutation is not selected at this stage, i.e. the initial sequence remains unchanged
    */
    Map_index_ener index_energies = do_step0(tabulist, Sref, vec_seq_residues, size_seq_PB, nbr_possibilities,
                                             seq_locked_pbe, montecarlo, threshold, threshold2);

    /* Print header information */
    string separator = "   ";
    cout << setfill(' ') << setw(5) << "Step" << separator
         << setfill(' ') << setw(5) << "Pos." << separator
         << setfill(' ') << setw(5) << "PBs" << separator
         << setfill(' ') << setw(Sref.get_all_blocks().size()+1) << "Sequence" << separator
         << setfill(' ') << setw(7) << "Energy" << separator
         << setfill(' ') << setw(8) << "TM-score" << separator
         << setfill(' ') << setw(8) << "GDT-TS" << separator
         << setfill(' ') << setw(6) << "RMSD" << endl;


    /* Step 1
       Select the PB to mutate in the initial sequence
       Main loop:
       - Sub loop on the tabu list to find the lowest pseudo-energy
       - The initial sequence is mutated by replacing the "b" PB at the "p" position corresponding to the lowest pseudo-energy
       - Update of the flag for the PB with the memory value of the tabu list
       - Calculate again all the energies (Step 0)
    */
    string iter;
    float new_energy, new_rmsd=0.00, new_tmscore=0.00, new_gdt=0.00, dU=0, rnd=0;
    float energy_ref=1000000, previous_energy=1000000; // arbitrarily high values

    /* Main loop */
    for (int i = 1; i <= steps; ++i) {
        /* Control the size of the index_energies map */
        if (index_energies.size() == 0){
            cerr << "* Program stopped: no more PB to test\n"
                 << "* Restart ISIDE with either:\n"
                 << "  (i) a smaller tabu list (-t arg)\n"
                 << "  (ii) a lower probability threshold (-r arg)\n"
                 << "  (iii) a larger number of PB to test (-n arg)\n"
                 << "  (iv) fewer locked positions (-l or -y arg)\n"
                 << endl;
            return 1;
        }

        /* Find lowest energy */
        new_energy = getMin(index_energies);

        if (montecarlo){
            if (i % decrease == 0) ck*=alpha;
            dU = previous_energy-new_energy;
            rnd = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        }

        if (!montecarlo || (montecarlo && rnd < exp (dU/ck))) {
            previous_energy = new_energy;
            vector<int> index_min_ener = getMinKey(index_energies);
            stringstream out, out2;
            string new_positions, new_blocks;
            PTools::Rigidbody new_struct = Sref;
    
            /* For each mutation */
            for (vector<int>::const_iterator j = index_min_ener.begin(); j != index_min_ener.end(); ++j) {
                /* Get the position and BP */
                int new_pos = tabulist[*j].get_pos();
                char new_block = tabulist[*j].get_block();
                /* Set -1 to the flag */
                tabulist.lock_mutation(new_pos, new_block);
    
                /* Mutate the structure */
                PTools::Rigidbody block_3D = sa.get_block(new_block);
                int residu_pos = new_pos + 2; /* always keep first 2 residus, do not start at 0 */
                new_struct = switch_PB(new_struct, residu_pos, block_3D);
    
                /* Cast new pos and new_block into string for printing later */
                out << residu_pos << " "; 
                new_blocks += string(1, new_block) + " ";
            }
    
            new_positions = out.str();
            new_struct.change_all_resnames(vec_seq_residues); /* Update residue names */
            new_struct.init_blocks_edges(); /* set ZZ at the begining and at the end */
    
            /* Save the new structure */
            out2 << setfill('0') << setw(5) << i; /* Cast iteration */
            iter = out2.str();
            string fileref = result_dirs[0] + "/" + iter + ".pdb";
            PTools::WritePDB(new_struct, fileref);

            /* Update RMSD and TM-score */
            if (aaseq.empty()){
                string command = string("./TMscore/TMscore") + " " + fileref + " " + pdbname;
                vector<string> rmsdtm = parseTM(&command[0u]);
                /* rmsdtm[0]: RMSD; rmsdtm[1]: TM-score; rmsdtm[2]: MaxSub-score;
                   rmsdtm[3]: GDT-TS; rmsdtm[4]: GDT-HA */
                new_rmsd = stof(rmsdtm[0]);
                new_tmscore = stof(rmsdtm[1]);
                new_gdt = stof(rmsdtm[3]);
            }

            if (montecarlo) remove(fileref.c_str()); /* Otherwise, too many files created */

            /* Print information */
            print_information(i, new_positions, new_blocks, new_struct, new_energy,
                              new_rmsd, printout, weblogo, new_tmscore, new_gdt);
        
            /* Update la liste Tabu with new energy */
            tabulist.update_ener(new_energy);
    
            /* Save structure with the best current energy */
            if (new_energy < energy_ref) {
                energy_ref = new_energy;
                PTools::WritePDB(new_struct, result_dirs[1] + "/" + iter + ".pdb");
            }
    
            /* Update variables */
            Sref = new_struct;
        }

        /* Perform the step 0 for the next mutation */
        index_energies = do_step0(tabulist, Sref, vec_seq_residues, size_seq_PB, nbr_possibilities,
                                  seq_locked_pbe, montecarlo, threshold, threshold2);
    } /* Main loop end */

    return 0;
}

