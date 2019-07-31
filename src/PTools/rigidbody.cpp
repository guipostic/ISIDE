/****************************************************************************
 *   Copyright (C) 2006-2008   Adrien Saladin                               *
 *   adrien.saladin@gmail.com                                               *
 *   Copyright (C) 2008   Pierre Poulain                                    *
 *   Copyright (C) 2008   Sebastien Fiorucci                                *
 *   Copyright (C) 2008   Chantal Prevost                                   *
 *   Copyright (C) 2008   Martin Zacharias                                  *
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


#include "rigidbody.h"
#include "atomselection.h"
#include "geometry.h"
#include "pdbio.h"

// ISIDE v1.0: Residues listing and compute functions
Residu::Residu(const std::string & type, uint _id){
    rtype=type;
    id=_id;
    pb='Z';
}


Residu::Residu(const std::string & type, uint _id, char _PB){
    rtype=type;
    id=_id;
    pb=_PB;
}


void Residu::set_atom(PTools::Atom atm){
    if(atm.GetType()=="CA") Ca=atm;
    if(atm.GetType()=="C") C=atm;
    if(atm.GetType()=="N") N=atm;
    if(atm.GetType()=="CB") Cb=atm;
}


double Residu::cb_angle(const PTools::Coord3D & coordAtm){
    PTools::Coord3D vect1;
    MakeVect(Ca.GetCoords(),Cb.GetCoords(),vect1);
    PTools::Coord3D vect2;
    MakeVect(Ca.GetCoords(),coordAtm,vect2);
    return(Angle(vect1,vect2));
}


void Residu::show(){
    std::cout<<id<< " " <<rtype<<" C : "<<C.GetCoords().x<<" "<<C.GetCoords().y<<" "<<C.GetCoords().z<<" Ca : "<<Ca.GetCoords().x<<" "<<Ca.GetCoords().y<<" "<<Ca.GetCoords().z<<" N : "<<N.GetCoords().x<<" "<<N.GetCoords().y<<" "<<N.GetCoords().z<< std::endl;
}


bool Residu::is_beta(){
    return(pb=='c' || pb=='d' || pb=='e');
}



namespace PTools{


Rigidbody::Rigidbody()
{
    ResetMatrix();
}



Rigidbody::Rigidbody(std::string filename)
{
    ReadPDB(filename,*this);
    ResetMatrix();
}

Rigidbody::Rigidbody(std::string filename, char block)
{
    ReadPDB(filename,*this, block);
    ResetMatrix();
}



Rigidbody::Rigidbody(const Rigidbody& model)
        : CoordsArray(model)
{
// ISIDE v1.0
// this copy constructor is needed because double[4][4] is not
// automatically copied with the default copy constructor
//TODO: check if always the case
    this->mAtomProp = model.mAtomProp;
    this-> _description = model._description;
}


void Rigidbody::AddAtom(const Atomproperty& at, Coord3D co)
{
    mAtomProp.push_back(at);
    AddCoord(co);
}


Atom Rigidbody::CopyAtom(uint i) const
{
    Atom at(mAtomProp[i],GetCoords(i));
    return at;
}

void Rigidbody::SetAtom(uint pos, const Atom& atom)
{
   //if (pos<0  || pos >= this->Size()) // ISIDE v1.0
   if (pos >= this->Size())
   {
      std::string message = "SetAtom: position ";
      message += pos;
      message += " is out of range";
      throw std::out_of_range(message);
   }
   Atomproperty atp(atom);
   Coord3D co(atom.GetCoords());
   SetAtomProperty(pos, atp);
   SetCoords(pos,co);
}



void Rigidbody::AddAtom(const Atom& at)
{
    Atomproperty atp(at);
    Coord3D co = at.GetCoords();
    AddAtom(atp,co);
}


Coord3D Rigidbody::FindCenter() const
{
    Coord3D center(0.0,0.0,0.0);
    for (uint i=0; i< this->Size() ; i++)
    {
        center =  center + GetCoords(i);
    }
    return ( (1.0/(double)this->Size())*center);
}


void Rigidbody::CenterToOrigin()
{
    Coord3D c = FindCenter();
    Translate(Coord3D()-c);
}


double Rigidbody::RadiusGyration()
{
    Coord3D c = this->FindCenter();
    double r=0.0;
    for (uint i=0; i< this->Size(); i++)
    {
        r += Norm2( c - this->GetCoords(i) );
    }

    double result = sqrt( r/ (double) this->Size() );
    return result;
}


double Rigidbody::Radius()
{
    Coord3D center = this->FindCenter();
    uint size = this->Size();
    double radius = 0.0;
    for (uint i=0; i < size; i++)
    {
        double rad=Norm(center - this->GetCoords(i));
        if (radius < rad) {radius=rad;}
    }
    return radius;
}


void Rigidbody::Translate(const Coord3D& tr)
{
    CoordsArray::Translate(tr);
}

void Rigidbody::AttractEulerRotate(double phi, double ssi, double rot)
{
   CoordsArray::AttractEulerRotate(phi, ssi, rot);
}


AtomSelection Rigidbody::SelectAllAtoms() const
{
    AtomSelection newsel;
    newsel.SetRigid(*this);
    for (uint i=0; i < Size(); i++)
    {
        newsel.AddAtomIndex(i);
    }

    return newsel;
}


AtomSelection Rigidbody::SelectAtomType(std::string atomtype)
{
    AtomSelection newsel;
    newsel.SetRigid(*this);

    for (uint i=0; i<Size(); i++)
    {
        if ( mAtomProp[i].GetType()==atomtype)
            newsel.AddAtomIndex(i);
    }

    return newsel;
}


AtomSelection Rigidbody::SelectResidType(std::string residtype)
{
    AtomSelection newsel;
    newsel.SetRigid(*this);

    for (uint i=0; i<Size(); i++)
    {
        if (mAtomProp[i].GetResidType()==residtype)
            newsel.AddAtomIndex(i);
    }
    return newsel;
}


AtomSelection Rigidbody::SelectChainId(std::string chainId) {
    AtomSelection newsel;
    newsel.SetRigid(*this);
    for (uint i=0; i<Size(); i++)
    {
        if (mAtomProp[i].GetChainId()==chainId)
            newsel.AddAtomIndex(i);
    }
    return newsel;
}


AtomSelection Rigidbody::SelectResRange(uint start, uint stop)
{
    AtomSelection newsel;
    newsel.SetRigid(*this);

    for (uint i=0; i < Size(); i++)
    {
        Atomproperty& atp ( mAtomProp[i] );
        if (atp.GetResidId() >=start && atp.GetResidId() <= stop) newsel.AddAtomIndex(i);
    }
    return newsel;
}


AtomSelection Rigidbody::CA() {
    return SelectAtomType("CA");
}


bool isBackbone(const std::string &  atomtype)
{

    const std::string bbtypes[] = {"N", "CA", "C", "O"};
    int const bbsize = sizeof(bbtypes)/sizeof(std::string);

    for (int i =0; i<bbsize; i++)
    {
        if (atomtype == bbtypes[i] ) return true;
    }

    return false;
}


AtomSelection Rigidbody::Backbone()
{
    AtomSelection newsel;
    newsel.SetRigid(*this);

    for (uint i=0; i<this->Size(); i++)
    {
        if (isBackbone(CopyAtom(i).GetType()) )
        {
            newsel.AddAtomIndex(i);
        }

    }
    return newsel;
}


/// operator +
Rigidbody Rigidbody::operator+(const Rigidbody& rig) {
    Rigidbody rigFinal(*this);
    for (uint i=0; i< rig.Size() ; i++) {
        rigFinal.AddCoord(rig.GetCoords(i));
        rigFinal.mAtomProp.push_back(rig.mAtomProp[i]);
    }
    return rigFinal;
}


void Rigidbody::ABrotate(const Coord3D& A, const Coord3D& B, double theta)
{
    PTools::ABrotate(A,B, *this, theta);
}


std::string Rigidbody::PrintPDB() const
{
    uint size=this->Size();

    std::string output;
    for (uint i=0; i < size ; i++)
    {
         Atom at(mAtomProp[i], this->GetCoords(i));
         output = output + at.ToPdbString();
    }
    return output;
}


void Rigidbody::ApplyMatrix(const Matrix& mat)
{

   double mat44[4][4];
   for(uint i=0; i<4;i++)
    for(uint j=0;j<4;j++)
      mat44[i][j]=mat(i,j);
   CoordsArray::MatrixMultiply(mat44);
}


/* ISIDE v1.0 */
void Rigidbody::Renumber(int start_pos){
    int previous_resid = start_pos - 1;
    int new_resid = previous_resid;
    for(uint i=0; i< mAtomProp.size(); i++){
        // For Atoms
        mAtomProp[i].SetAtomId(i+1);

        // For residues
        if ( previous_resid != (int)mAtomProp[i].GetResidId() ){
            new_resid = new_resid + 1;
            previous_resid = mAtomProp[i].GetResidId();
        }
        mAtomProp[i].SetResidId(new_resid);
    }
}


void Rigidbody::Renumber_res(){
    int nbr_atom_per_residue = 4;
    int resid = 0;
    for(uint i=0; i< mAtomProp.size(); i++){
        // For residues
        if (i%nbr_atom_per_residue == 0){
            resid = resid + 1;
        }
        mAtomProp[i].SetResidId(resid);
    }
}


/* ISIDE v1.0 */
double Rigidbody::compute_energies(const Ppot &d, Solvation sol, Bsheet beta){
    float angleBeta=0;
    float distance=0;
    double eppot=0,esol=0,ebsheet=0;
    int nbr_atoms=0;
    size_t temp_resid=0;
    std::string res_name1, res_name2, atom_name1, atom_name2,key;
    Atom atom1,atom2;
    int resid1=0,resid2=0;
    Coord3D vect1;
    Coord3D vect2;
    std::vector<Residu> mResidu;
    this->list_res(mResidu);    

    for (uint i=0;i<Size();i++){
        //Get properties of atom i
        atom1 = GetAtom(i);
        res_name1 = atom1.GetResidType();
        atom_name1 = atom1.GetType();
        if (atom_name1=="CA")temp_resid=atom1.GetResidId();
        for (uint j=0;j<Size();j++){
            //Get properties of atom j
            atom2 = GetAtom(j);
            res_name2 = atom2.GetResidType();
            atom_name2 = atom2.GetType();
            resid1=(int)atom1.GetResidId();
            resid2=(int)atom2.GetResidId();
            
            // Considers only residues separated by 4 posiitons
            //if(atom1.GetResidId() > atom2.GetResidId()+4 || atom1.GetResidId() < atom2.GetResidId()-4)
            if (resid1-resid2>4 || resid1-resid2<-4){
                if (atom1.GetAtomId() == 1 && atom_name1.substr(0,1) == "N") { continue; }
                //if(atom2.GetAtomId() == 1 && atom_name2.substr(0,1) == "N") { continue; }
                
                key=res_name1+atom_name1+res_name2+atom_name2;

                // Exception if pairwise atomic interaction not described in the Ppot.par file
                distance=Dist(atom1,atom2);
                try {
                    Id a = d[key];
                    if (distance<15){
                        eppot+=d.getEnergy(a,distance);
                    }

                } catch (std::exception &e) {
                    //cout << "key : " << key << endl;
                }
                if (atom_name1 == "CA"){                
                    if (atom_name2 == "CA"){
                        //======== Hydropobicity part
                        if (distance<12 && mResidu[temp_resid-1].cb_angle(atom2.GetCoords())<1.570796)nbr_atoms+=1;
                        //======== Beta sheets part
                        // The 2 residues must correspond to PBs c, d, or e
                        if (mResidu[temp_resid-1].is_beta() && mResidu[atom2.GetResidId()-1].is_beta()){
                            // Do not applies to the residues at the extremities (look at i-1 and i+1)
                            if (temp_resid>1 && temp_resid<Size() && atom2.GetResidId()>1 && atom2.GetResidId()<Size()){
                                MakeVect(mResidu[temp_resid-2].Ca.GetCoords(),mResidu[temp_resid].Ca.GetCoords(),vect1);
                                MakeVect(mResidu[atom2.GetResidId()-2].Ca.GetCoords(),mResidu[atom2.GetResidId()].Ca.GetCoords(),vect2);
                                
                                angleBeta=Angle(vect1,vect2)*(180.0/3.14159265359);
                                
                                if (isnanf(angleBeta))
                                    std::cout<<vect1.x<<" "<<vect1.y<<" "<<vect1.z<<" "<<vect2.x<<" "<<vect2.y<<" "<<vect2.z<<std::endl;

                                ebsheet+=beta.getEnergy(distance,angleBeta);
                            }
                        }
                    }
                }
            }
        }
        if (atom_name1 == "CA"){
            esol+=sol.getEnergy(res_name1,nbr_atoms);
            nbr_atoms=0;
        }
    }

    return (esol+eppot+ebsheet);
}


std::vector<std::string> Rigidbody::get_all_resnames() {
    std::vector<std::string> list_resnames;
    for(size_t i=0; i < mAtomProp.size(); i++)
        if(i == 0 || mAtomProp[i].GetResidId() != mAtomProp[i-1].GetResidId())
            list_resnames.push_back(mAtomProp[i].GetResidType());

    return list_resnames;
}


void Rigidbody::change_all_resnames(const std::vector<std::string> & resnames) {
    int res_pos = 0;
    for(size_t i=0; i < mAtomProp.size(); i++) {
        // Get the change of residue
        if (i != 0 && mAtomProp[i].GetResidId() != mAtomProp[i-1].GetResidId()) { res_pos++; }
        mAtomProp[i].SetResidType(resnames[res_pos]);
    }
}


void Rigidbody::change_resname(const std::string &resname, int pos) {
    mAtomProp[pos].SetResidType(resname);
}


// Merge 2 blocks after block2 superposed to block1
Rigidbody Merge(const Rigidbody &block1, const Rigidbody &block2, uint start_block1, uint end_block1, uint start_block2, uint end_block2)
{
    if (end_block1 > block1.Size() || start_block1 > end_block1 || end_block2 > block2.Size() || start_block2 > end_block2 )
        throw ("Not a good index");

    //Need tmp because of the "const"
    Rigidbody tmp1 = block1;
    Rigidbody tmp2 = block2;
    Rigidbody result;

    //PTools::AtomSelection selected_atom_block1(tmp1.SelectResRange(1,(tmp1.Size()/4)-2));
    AtomSelection selected_atom_block1(tmp1.SelectResRange(start_block1,end_block1));
    Rigidbody r1 = selected_atom_block1.CreateRigid();

    //PTools::AtomSelection selected_atom_block2(tmp2.SelectResRange(3,5));
    AtomSelection selected_atom_block2(tmp2.SelectResRange(start_block2,end_block2));
    Rigidbody r2 = selected_atom_block2.CreateRigid();

    r2.Renumber();
    /* PTools::WritePDB(r2,"r2.pdb"); */

    result = r1 + r2;

    return result;
}


std::string Rigidbody::get_all_blocks(bool print_pos) const {

    std::string blocks = "";
    for(uint i=0; i < mAtomProp.size() ; i++) {
        if(i == 0 || mAtomProp[i].GetResidId() != mAtomProp[i-1].GetResidId()) {
            std::stringstream out2;
            if(print_pos) {
                out2 <<  mAtomProp[i].GetBlock() << "(" << mAtomProp[i].GetBlockInPos() << ")"; 
            } else {
                out2 <<  mAtomProp[i].GetBlock();
            }
            std::string s = out2.str();
            blocks+= s;
        }
    }
    return blocks;
}


char Rigidbody::get_block(int pos_res) const {
    for(uint i = 0; i < mAtomProp.size() ; i++)
        if ( pos_res == (int)mAtomProp[i].GetResidId() )
            return mAtomProp[i].GetBlock();

    return 'Z';
}


void Rigidbody::set_block_residue(char block, int pos_res) {
    for(uint i = 0; i < mAtomProp.size(); i++)
        if ( pos_res == (int)mAtomProp[i].GetResidId() )
            mAtomProp[i].SetBlock(block);
}


/* ISIDE v1.0 */
void Rigidbody::set_all_blocks(const std::string &blocks) {
    size_t res_pos = 0;
    for(uint i = 0; i < mAtomProp.size() ; i++) {
        // Get the change of residu
        if(i != 0 && mAtomProp[i].GetResidId() != mAtomProp[i-1].GetResidId()) {
            res_pos++;
        }
        if(res_pos != 0 && res_pos != 1 && res_pos != blocks.size()+4-1 && res_pos != blocks.size()+4-2) {
            mAtomProp[i].SetBlock(blocks[res_pos-2]);
        }
    }
}


/* ISIDE v1.0 */
void Rigidbody::set_all_blocks(char block) {
    for(uint i = 0; i < mAtomProp.size() ; i++)
        mAtomProp[i].SetBlock(block);
}


/* ISIDE v1.0 */
void Rigidbody::init_blocks_edges() {
    //TODO: a améliorer pour eviter le boucle sur tous les atomes
    size_t res_pos = 0;
    for(uint i = 0; i < mAtomProp.size() ; i++) {
        // Get the change of residu
        if(i != 0 && mAtomProp[i].GetResidId() != mAtomProp[i-1].GetResidId())
            res_pos++;

        if(res_pos == 0 || res_pos == 1 || res_pos == (this->Size_residus()-2) || res_pos == this->Size_residus()-1 )
            mAtomProp[i].SetBlock('Z');
    }  

}


/* ISIDE v1.0 */
void Rigidbody::set_inpos_block_per_residue(int inpos, int pos_res) {
    for(uint i = 0; i < mAtomProp.size() ; i++)
        if ( pos_res == (int)mAtomProp[i].GetResidId() )
            mAtomProp[i].SetBlockInPos(inpos);
}


std::string Rigidbody::get_inpos_blocks() const {
    std::string intern_pos = "";
    for(uint i=0; i < mAtomProp.size() ; i++) {
        if(i == 0 || mAtomProp[i].GetResidId() != mAtomProp[i-1].GetResidId()) {
            std::stringstream out2;
            out2 <<  mAtomProp[i].GetBlockInPos();
            intern_pos += out2.str();
        }
    }
    return intern_pos;  
}


void Rigidbody::list_res(std::vector<Residu>& rlist ){
    std::string string_pb_seq=this->get_all_blocks(false);
    const char* t_pb_seq=string_pb_seq.c_str();
    //std::cout<<"plop:"<<pb_seq<<"\n";
    int res_pos = 0;
    rlist.push_back(Residu(mAtomProp[0].GetResidType(),mAtomProp[0].GetResidId(),t_pb_seq[0]));
    for(uint i = 0; i < mAtomProp.size() ; i++) {
        //Get the change of residu
        if(i != 0 && mAtomProp[i].GetResidId() != mAtomProp[i-1].GetResidId()) {
            res_pos++;
            rlist.push_back(Residu(mAtomProp[i].GetResidType(),mAtomProp[i].GetResidId(),t_pb_seq[res_pos]));
        }
        rlist[res_pos].set_atom(GetAtom(i));
    }
}

} //namespace PTools
