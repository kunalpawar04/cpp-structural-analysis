#include <bits/stdc++.h>
using namespace std;

class parameters{

public:

    float sigma_x,sigma_y,sigma_z,Txy,Tyz,Txz,nx,ny,nz,i1,i2,i3;

    void input(int ip){

        cout<<"\n INPUT:\n";

        cout<<"\n Normal Stress x = ";
        cin>>sigma_x;

        cout<<"\n Normal Stress y = ";
        cin>>sigma_y;

        cout<<"\n Normal Stress z = ";
        cin>>sigma_z;

        cout<<"\n Shear Stress XY = ";
        cin>>Txy;

        cout<<"\n Shear Stress YZ = ";
        cin>>Tyz;

        cout<<"\n Shear Stress XZ = ";
        cin>>Txz;

        if(ip==1){

            cout<<"\n {If the direction cosine is not known give -1 as input}\n";

            while(true){

                cout<<"\n Enter direction Cosines: ";
                cout<<"\n nx: ";
                cin>>nx;
                cout<<"\n ny: ";
                cin>>ny;
                cout<<"\n nz: ";
                cin>>nz;
                int count_of_one=0;
                if(nx==-1)count_of_one++;
                if(ny==-1)count_of_one++;
                if(nz==-1)count_of_one++;
                if(count_of_one>1){
                    cout<<"\n You need to enter at least two direction cosines";
                    continue;
                }
                if(count_of_one==0){
                    if(nx*nx+ny*ny+nz*nz>1){
                        cout<<"\n The value nx*nx+ny*ny+nz*nz cannot be greater than 1";
                        continue;
                    }
                    break;
                }
                if(nx==-1){
                    if(ny*ny+nz*nz>1){
                        cout<<"\n The value ny*ny+nz*nz cannot be greater than 1";
                        continue;
                    }
                    nx=sqrt(1-ny*ny-nz*nz);
                    break;
                }
                if(ny==-1){
                    if(nx*nx+nz*nz>1){
                        cout<<"\n The value nx*nx+nz*nz cannot be greater than 1";
                        continue;
                    }
                    ny=sqrt(1-nx*nx-nz*nz);
                    break;
                }
                if(nz==-1){
                    if(ny*ny+nx*nx>1){
                        cout<<"\n The value ny*ny+nx*nx cannot be greater than 1";
                        continue;
                    }
                    nz=sqrt(1-nx*nx-ny*ny);
                    break;
                }
            }
        }
        cout<<"\n \n OUTPUT: \n ";
    }
    void invarients(){

        i1 = -(sigma_x + sigma_y + sigma_z);
        i2 = sigma_x*sigma_y + sigma_x*sigma_z + sigma_y*sigma_z - Txy*Txy - Tyz*Tyz - Txz*Txz;
        i3 = -(sigma_x*sigma_y*sigma_z + 2*Txy*Tyz*Txz - sigma_x*Tyz*Tyz - sigma_y*Txz*Txz - sigma_z*Txy*Txy);

    }

};
class stress_on_arbitrary_plane: public parameters{

    float Tres,sigma,Tx,Ty,Tz,T;
public:

    void calculate(){

        Tx = nx*sigma_x + ny*Txy + nz*Txz;
        Ty = nx*Txy + ny*sigma_y + nz*Tyz;
        Tz = ny*Txz + ny*Tyz + nz*sigma_z;

    }
    float resultant(){

        Tres = sqrt(Tx*Tx + Ty*Ty + Tz*Tz);
        return Tres;

    }
    float normal_stress(){

        sigma = nx*Tx+ny*Ty+nz*Tz;
        return sigma;

    }
    float shear_stress(){

        T = sqrt(Tres*Tres - sigma*sigma);
        return T;

    }
};
class stress_on_principal_plane: public parameters{

    float a,b,phi,g,A,B,C,nx,ny,nz;
    vector<float>sigma;

public:
    void calculate(){

        invarients();
        a = (3*i2 - i1*i1)/3.0;
        b = (2*i1*i1*i1 - 9*i1*i2 + 27*i3)/27.0;
        phi = acos(-b/(2.0*sqrt(-(a*a*a)/27.0)));
        g = 2.0*sqrt(-a/3.0);

    }
    void stresses(){

        sigma.resize(3);
        sigma[0] = g*cos(phi/3.0) - i1/3.0;
        sigma[1] = g*cos((phi + 2.0*acos(0.0)*2.0)/3.0) - i1/3.0;
        sigma[2] = g*cos((phi + 2.0*acos(0.0)*4.0)/3.0) - i1/3.0;


        for(int i = 0 ; i < 3 ; i++)cout<<" Principal Stress "<<i+1<<" : "<<sigma[i]<<"\n";
        sort(sigma.rbegin(),sigma.rend());

    }
    void directions(){

        sigma_x -= sigma[0];
        sigma_y -= sigma[0];
        sigma_z -= sigma[0];

        A = sigma_y*sigma_z - Tyz*Tyz;
        B = Txy*sigma_z - Txz*Tyz;
        C = Txy*Tyz - Txz*sigma_y;

        float uT = sqrt( A*A + B*B + C*C);

        nx = abs(A/uT);
        ny = abs(B/uT);
        nz = abs(C/uT);

        cout<<" "<<nx<<"\n "<<ny<<"\n "<<nz<<"\n ";

    }
};
class stress_on_octahedral_plane:public parameters{

    float sigma_octahedral,tau_octahedral;

public:

    void calculate(){

        invarients();

        sigma_octahedral=abs(i1)/3.0;
        tau_octahedral=sqrt(2*i1*i1-6*abs(i2))/3.0;

    }
    void get_stress(){

        cout<<"\n Normal Stress on Octahedral plane: "<<sigma_octahedral<<'\n';
        cout<<"\n Shear stress on Octahedral plane: "<<tau_octahedral<<'\n';

    }
};
int main()
{

    parameters par;

    stress_on_arbitrary_plane plane_stress;

    stress_on_principal_plane principal_stress;

    stress_on_octahedral_plane octaherdal_stress;

    while(true){

        cout<<"\n Enter 1 to find \"Stress on an Arbitrary plane\"\n";

        cout<<" Enter 2 to find \"Stress on Principal Plane\"\n";

        cout<<" Enter 3 to find \"Stress on Octahedral Plane\"\n";

        cout<<" Enter 0 to terminate program\n ";

        int input;
        cin>>input;

        if ( input==1 ){

            plane_stress.input(input);

            plane_stress.calculate();

            cout<<" Resultant Stress = "<<plane_stress.resultant()<<"\n";

            cout<<" Normal Stress = "<<plane_stress.normal_stress()<<"\n";

            cout<<" Shear Stress = "<<plane_stress.shear_stress()<<"\n";

        }
        else if ( input==2 ){

            principal_stress.input(input);

            principal_stress.calculate();

            cout<<" Normal Stresses on Principal Plane are: \n";
            principal_stress.stresses();

            cout<<" Directions of Principal Plane are: \n";
            principal_stress.directions();

        }
        else if ( input ==3 ){

            octaherdal_stress.input(input);

            octaherdal_stress.calculate();

            octaherdal_stress.get_stress();
        }
        else break;
    }
    return 0;
}
