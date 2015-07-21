#include "OpenChannel3D.h"
#include "lattice_vars.h"
#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
#include "workArounds.h"

// for debugging
#include <iostream>

using namespace std;

OpenChannel3D::OpenChannel3D(const int rk, const int sz, const string input_file):
rank(rk), size(sz)
{
    tag_d = 666; tag_u = 999;
    read_input_file(input_file);
    initialize_lattice_data();
    initialize_local_partition_variables();
    initialize_mpi_buffers();
    vtk_ts = 0;
}

OpenChannel3D::~OpenChannel3D(){
    delete [] fEven;
    delete [] fOdd;
    delete [] inl;
    delete [] onl;
    delete [] snl;
    delete [] u_bc;
    
    delete [] ghost_in_m;
    delete [] ghost_out_m;
    delete [] ghost_in_p;
    delete [] ghost_out_p;
    
    delete [] rho_l;
    delete [] ux_l;
    delete [] uy_l;
    delete [] uz_l;
    
}

void OpenChannel3D::write_data(MPI_Comm comm, bool isEven){
    
    string densityFileStub("density");
    string ux_FileStub("ux");
    string uy_FileStub("uy");
    string uz_FileStub("uz");
    stringstream ts_ind;
    string ts_ind_str;
    string fileSuffix(".b_dat");
    string density_fn, ux_fn, uy_fn, uz_fn;
    MPI_File fh_rho, fh_ux, fh_uy, fh_uz;
    MPI_Status mpi_s1, mpi_s2, mpi_s3,mpi_s4;
    
    //	// create temporary data buffers
    //	float * rho_l = new float[numEntries];
    //	float * ux_l = new float[numEntries];
    //	float * uy_l = new float[numEntries];
    //	float * uz_l = new float[numEntries];
    
    float * fOut;
    if (isEven){
        fOut = fEven;
    }else{
        fOut = fOdd;
    }
    
    int nnodes = this->nnodes;
    int numSpd = this->numSpd;
    int Nx = this->Nx;
    int Ny = this->Ny;
    const float * ex = this->ex;
    const float * ey = this->ey;
    const float * ez = this->ez;
    float * ux_l = this->ux_l;
    float * uy_l = this->uy_l;
    float * uz_l = this->uz_l;
    float * rho_l = this->rho_l;
    int totalSlices = this->totalSlices;
    int * snl = this->snl;
    int numMySlices = this->numMySlices;
    int numEntries = Nx*Ny*numMySlices;
    dummyUse(nnodes,numEntries);
    
    #pragma omp parallel for collapse(3)
    #pragma acc parallel loop collapse(3) \
        present(fOut[0:nnodes*numSpd], snl[0:nnodes]) \
        copyout(ux_l[0:numEntries],uy_l[0:numEntries],uz_l[0:numEntries],rho_l[0:numEntries]) \
        copyin(ex[0:numSpd],ey[0:numSpd],ez[0:numSpd])
    for(int z = HALO;z<(totalSlices-HALO);z++){
        for(int y = 0;y<Ny;y++){
            for(int x = 0;x<Nx;x++){
                int tid_l, tid_g;
                float tmp_rho, tmp_ux, tmp_uy, tmp_uz;
                tid_l = x+y*Nx+(z-HALO)*Nx*Ny;
                tmp_rho = 0; tid_g = x+y*Nx+z*Nx*Ny;
                tmp_ux = 0; tmp_uy = 0; tmp_uz = 0;
                for(int spd=0;spd<numSpd;spd++){
                    float f = fOut[getIdx(nnodes,numSpd,tid_g,spd)];
                    tmp_rho+=f;
                    tmp_ux+=ex[spd]*f;
                    tmp_uy+=ey[spd]*f;
                    tmp_uz+=ez[spd]*f;
                }
                rho_l[tid_l]=tmp_rho;
                if(snl[tid_g]==1){
                    ux_l[tid_l]=0.; uy_l[tid_l]=0.; uz_l[tid_l]=0.;
                }else{
                    ux_l[tid_l]=tmp_ux*(1./tmp_rho);
                    uy_l[tid_l]=tmp_uy*(1./tmp_rho);
                    uz_l[tid_l]=tmp_uz*(1./tmp_rho);
                }
            }
        }
    }
    
    // generate file names
    ts_ind << vtk_ts;
    density_fn = densityFileStub+ts_ind.str()+fileSuffix;
    ux_fn = ux_FileStub+ts_ind.str()+fileSuffix;
    uy_fn = uy_FileStub+ts_ind.str()+fileSuffix;
    uz_fn = uz_FileStub+ts_ind.str()+fileSuffix;
    
    // open MPI file for parallel IO
    MPI_File_open(comm,(char*)density_fn.c_str(),
    MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh_rho);
    
    MPI_File_open(comm,(char*)ux_fn.c_str(),
    MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh_ux);
    
    MPI_File_open(comm,(char*)uy_fn.c_str(),
    MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh_uy);
    
    MPI_File_open(comm,(char*)uz_fn.c_str(),
    MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh_uz);
    
    //write your chunk of data
    MPI_File_write_at(fh_rho,offset,rho_l,numEntries,MPI_FLOAT,&mpi_s1);
    MPI_File_write_at(fh_ux,offset,ux_l,numEntries,MPI_FLOAT,&mpi_s2);
    MPI_File_write_at(fh_uy,offset,uy_l,numEntries,MPI_FLOAT,&mpi_s3);
    MPI_File_write_at(fh_uz,offset,uz_l,numEntries,MPI_FLOAT,&mpi_s4);
    
    //close the files
    MPI_File_close(&fh_rho);
    MPI_File_close(&fh_ux);
    MPI_File_close(&fh_uy);
    MPI_File_close(&fh_uz);
    
    vtk_ts++; // increment the dump counter...
}

void OpenChannel3D::D3Q15_process_slices(bool isEven, const int firstSlice, const int lastSlice){
    
    // this monstrosity needs to be change into something more simple and clear.
    // it performs on the GPU but is rather unmaintainable.
    
    const float * fIn;
    float * fOut;
    
    if(isEven){
        fIn = fEven; fOut = fOdd;
    }else{
        fIn = fOdd; fOut = fEven;
    }
    
    // local copies of class data members needed for acc compiler
    int* inl = this->inl;
    int* onl = this->onl;
    int* snl = this->snl;
    float* u_bc = this->u_bc;
    int Ny = this->Ny;
    int Nx = this->Nx;
    float omega = this->omega;
    int nnodes = this->nnodes;
    
    dummyUse(nnodes);
    
    
    //Nz=lastSlice-firstSlice;
    const int numSpd=15;
    #pragma omp parallel for collapse(3)
    #pragma acc parallel loop collapse(3) \
        present(fIn[0:nnodes*numSpd]) \
        present(fOut[0:nnodes*numSpd]) \
        present(inl[0:nnodes], onl[0:nnodes], snl[0:nnodes], u_bc[0:nnodes])
    for(int Z=firstSlice;Z<lastSlice;Z++){
        for(int Y=0;Y<Ny;Y++){
            for(int X=0;X<Nx;X++){
                float f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14;
                float cu,rho,ux,uy,uz,fEq,dz;
                int X_t,Y_t,Z_t,tid_t,tid;
                
                tid=X+Y*Nx+Z*Nx*Ny;
                
                //load the data into registers
                // f0=fIn[tid]; f1=fIn[Nx*Ny*Nz+tid];
                // f2=fIn[2*Nx*Ny*Nz+tid]; f3=fIn[3*Nx*Ny*Nz+tid];
                // f4=fIn[4*Nx*Ny*Nz+tid]; f5=fIn[5*Nx*Ny*Nz+tid];
                // f6=fIn[6*Nx*Ny*Nz+tid]; f7=fIn[7*Nx*Ny*Nz+tid];
                // f8=fIn[8*Nx*Ny*Nz+tid]; f9=fIn[9*Nx*Ny*Nz+tid];
                // f10=fIn[10*Nx*Ny*Nz+tid]; f11=fIn[11*Nx*Ny*Nz+tid];
                // f12=fIn[12*Nx*Ny*Nz+tid]; f13=fIn[13*Nx*Ny*Nz+tid];
                // f14=fIn[14*Nx*Ny*Nz+tid];
                f0=fIn[getIdx(nnodes, numSpd, tid,0)]; f1=fIn[getIdx(nnodes, numSpd, tid,1)];
                f2=fIn[getIdx(nnodes, numSpd, tid,2)]; f3=fIn[getIdx(nnodes, numSpd, tid,3)];
                f4=fIn[getIdx(nnodes, numSpd, tid,4)]; f5=fIn[getIdx(nnodes, numSpd, tid,5)];
                f6=fIn[getIdx(nnodes, numSpd, tid,6)]; f7=fIn[getIdx(nnodes, numSpd, tid,7)];
                f8=fIn[getIdx(nnodes, numSpd, tid,8)]; f9=fIn[getIdx(nnodes, numSpd, tid,9)];
                f10=fIn[getIdx(nnodes, numSpd, tid,10)]; f11=fIn[getIdx(nnodes, numSpd, tid,11)];
                f12=fIn[getIdx(nnodes, numSpd, tid,12)]; f13=fIn[getIdx(nnodes, numSpd, tid,13)];
                f14=fIn[getIdx(nnodes, numSpd, tid,14)];
                
                //compute density
                rho = f0+f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f11+f12+f13+f14;
                ux=f1-f2+f7-f8+f9-f10+f11-f12+f13-f14; ux/=rho;
                uy=f3-f4+f7+f8-f9-f10+f11+f12-f13-f14; uy/=rho;
                uz=f5-f6+f7+f8+f9+f10-f11-f12-f13-f14; uz/=rho;
                
                //if it's on the inl or onl, update
                
                if((inl[tid]==1)||(onl[tid]==1)){
                    
                    dz=u_bc[tid]-uz;
                    //speed 1 ex=1 ey=ez=0. w=1./9.
                    cu=3.F*(1.F)*(-ux);
                    f1+=(1.F/9.F)*rho*cu;
                    
                    //speed 2 ex=-1 ey=ez=0. w=1./9.
                    cu=3.F*(-1.F)*(-ux);
                    f2+=(1.F/9.F)*rho*cu;
                    
                    //speed 3 ey=1; ex=ez=0; w=1./9.
                    cu=3.F*(1.F)*(-uy);
                    f3+=(1.F/9.F)*rho*cu;
                    
                    //speed 4 ey=-1; ex=ez=0; w=1./9.
                    cu=3.F*(-1.F)*(-uy);
                    f4+=(1.F/9.F)*rho*cu;
                    
                    //speed 5 ex=ey=0; ez=1; w=1./9.
                    cu=3.F*(1.F)*(dz);
                    f5+=(1.F/9.F)*rho*cu;
                    
                    //speed 6 ex=ey=0; ez=-1; w=1./9.
                    cu=3.F*(-1.F)*(dz);
                    f6+=(1.F/9.F)*rho*cu;
                    
                    //speed 7 ex=ey=ez=1; w=1./72.
                    cu=3.F*((1.F)*-ux+(1.F)*(-uy)+(1.F)*dz);
                    f7+=(1.F/72.F)*rho*cu;
                    
                    //speed 8 ex=-1 ey=ez=1; w=1./72.
                    cu=3.F*((-1.F)*-ux+(1.F)*(-uy)+(1.F)*dz);
                    f8+=(1.F/72.F)*rho*cu;
                    
                    //speed 9 ex=1 ey=-1 ez=1
                    cu=3.0F*((1.F)*-ux+(-1.F)*(-uy)+(1.F)*dz);
                    f9+=(1.F/72.F)*rho*cu;
                    
                    //speed 10 ex=-1 ey=-1 ez=1
                    cu=3.0F*((-1.F)*-ux+(-1.F)*(-uy)+(1.F)*dz);
                    f10+=(1.F/72.F)*rho*cu;
                    
                    //speed 11 ex=1 ey=1 ez=-1
                    cu=3.0F*((1.F)*-ux +(1.F)*(-uy)+(-1.F)*dz);
                    f11+=(1.F/72.F)*rho*cu;
                    
                    //speed 12 ex=-1 ey=1 ez=-1
                    cu=3.0F*((-1.F)*-ux+(1.F)*(-uy)+(-1.F)*dz);
                    f12+=(1.F/72.F)*rho*cu;
                    
                    //speed 13 ex=1 ey=-1 ez=-1 w=1./72.
                    cu=3.0F*((1.F)*-ux+(-1.F)*(-uy)+(-1.F)*dz);
                    f13+=(1.F/72.F)*rho*cu;
                    
                    //speed 14 ex=ey=ez=-1 w=1./72.
                    cu=3.0F*((-1.F)*-ux + (-1.F)*(-uy) +(-1.F)*dz);
                    f14+=(1.F/72.F)*rho*cu;
                    
                    ux=0.; uy=0.; uz=u_bc[tid];
                }
                
                if(snl[tid]==1){
                    // 1--2
                    cu=f2; f2=f1; f1=cu;
                    //3--4
                    cu=f4; f4=f3; f3=cu;
                    //5--6
                    cu=f6; f6=f5; f5=cu;
                    //7--14
                    cu=f14; f14=f7; f7=cu;
                    //8--13
                    cu=f13; f13=f8; f8=cu;
                    //9--12
                    cu=f12; f12=f9; f9=cu;
                    //10--11
                    cu=f11; f11=f10; f10=cu;
                    
                    
                }else{
                    fEq=rho*(2.F/9.F)*(1.F-1.5F*(ux*ux+uy*uy+uz*uz));
                    f0=f0-omega*(f0-fEq);
                    
                    //speed 1 ex=1 ey=ez=0 w=1./9.
                    cu=3.F*(1.F*ux);
                    fEq=rho*(1.F/9.F)*(1.F+cu+0.5F*(cu*cu)-
                    1.5F*(ux*ux+uy*uy+uz*uz));
                    f1=f1-omega*(f1-fEq);
                    
                    //speed 2 ex=-1 ey=ez=0 w=1./9.
                    cu=3.F*((-1.F)*ux);
                    fEq=rho*(1.F/9.F)*(1.F+cu+0.5F*(cu*cu)-
                    1.5F*(ux*ux+uy*uy+uz*uz));
                    f2=f2-omega*(f2-fEq);
                    
                    //speed 3 ex=0 ey=1 ez=0 w=1./9.
                    cu=3.F*(1.F*uy);
                    fEq=rho*(1.F/9.F)*(1.F+cu+0.5F*(cu*cu)-
                    1.5F*(ux*ux+uy*uy+uz*uz));
                    f3=f3-omega*(f3-fEq);
                    
                    //speed 4 ex=0 ey=-1 ez=0 w=1./9.
                    cu=3.F*(-1.F*uy);
                    fEq=rho*(1.F/9.F)*(1.F+cu+0.5F*(cu*cu)-
                    1.5F*(ux*ux+uy*uy+uz*uz));
                    f4=f4-omega*(f4-fEq);
                    
                    //speed 5 ex=ey=0 ez=1 w=1./9.
                    cu=3.F*(1.F*uz);
                    fEq=rho*(1.F/9.F)*(1.F+cu+0.5F*(cu*cu)-
                    1.5F*(ux*ux+uy*uy+uz*uz));
                    f5=f5-omega*(f5-fEq);
                    
                    //speed 6 ex=ey=0 ez=-1 w=1./9.
                    cu=3.F*(-1.F*uz);
                    fEq=rho*(1.F/9.F)*(1.F+cu+0.5F*(cu*cu)-
                    1.5F*(ux*ux+uy*uy+uz*uz));
                    f6=f6-omega*(f6-fEq);
                    
                    //speed 7 ex=ey=ez=1 w=1./72.
                    cu=3.F*(ux+uy+uz);
                    fEq=rho*(1.F/72.F)*(1.F+cu+0.5F*(cu*cu)-
                    1.5F*(ux*ux+uy*uy+uz*uz));
                    f7=f7-omega*(f7-fEq);
                    
                    //speed 8 ex=-1 ey=ez=1 w=1./72.
                    cu=3.F*(-ux+uy+uz);
                    fEq=rho*(1.F/72.F)*(1.F+cu+0.5F*(cu*cu)-
                    1.5F*(ux*ux+uy*uy+uz*uz));
                    f8=f8-omega*(f8-fEq);
                    
                    //speed 9 ex=1 ey=-1 ez=1 w=1./72.
                    cu=3.F*(ux-uy+uz);
                    fEq=rho*(1.F/72.F)*(1.F+cu+0.5F*(cu*cu)-
                    1.5F*(ux*ux+uy*uy+uz*uz));
                    f9=f9-omega*(f9-fEq);
                    
                    //speed 10 ex=-1 ey=-1 ez=1 w=1/72
                    cu=3.F*(-ux-uy+uz);
                    fEq=rho*(1.F/72.F)*(1.F+cu+0.5F*(cu*cu)-
                    1.5F*(ux*ux+uy*uy+uz*uz));
                    f10=f10-omega*(f10-fEq);
                    
                    //speed 11 ex=1 ey=1 ez=-1 w=1/72
                    cu=3.F*(ux+uy-uz);
                    fEq=rho*(1.F/72.F)*(1.F+cu+0.5F*(cu*cu)-
                    1.5F*(ux*ux+uy*uy+uz*uz));
                    f11=f11-omega*(f11-fEq);
                    
                    //speed 12 ex=-1 ey=1 ez=-1 w=1/72
                    cu=3.F*(-ux+uy-uz);
                    fEq=rho*(1.F/72.F)*(1.F+cu+0.5F*(cu*cu)-
                    1.5F*(ux*ux+uy*uy+uz*uz));
                    f12=f12-omega*(f12-fEq);
                    
                    //speed 13 ex=1 ey=ez=-1 w=1/72
                    cu=3.F*(ux-uy-uz);
                    fEq=rho*(1.F/72.F)*(1.F+cu+0.5F*(cu*cu)-
                    1.5F*(ux*ux+uy*uy+uz*uz));
                    f13=f13-omega*(f13-fEq);
                    
                    //speed 14 ex=ey=ez=-1 w=1/72
                    cu=3.F*(-ux-uy-uz);
                    fEq=rho*(1.F/72.F)*(1.F+cu+0.5F*(cu*cu)-
                    1.5F*(ux*ux+uy*uy+uz*uz));
                    f14=f14-omega*(f14-fEq);
                    
                    
                    
                }
                
                //speed 0 ex=ey=ez=0
                //fOut[tid]=f0;
                fOut[getIdx(nnodes, numSpd, tid,0)]=f0;
                
                //speed 1 ex=1 ey=ez=0
                X_t=X+1; Y_t=Y; Z_t=Z;
                if(X_t==Nx) X_t=0;
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[Nx*Ny*Nz+tid_t]=f1;
                fOut[getIdx(nnodes, numSpd, tid_t,1)]=f1;
                
                //speed 2 ex=-1 ey=ez=0;
                X_t=X-1; Y_t=Y; Z_t=Z;
                if(X_t<0) X_t=(Nx-1);
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[2*Nx*Ny*Nz+tid_t]=f2;
                fOut[getIdx(nnodes, numSpd, tid_t,2)]=f2;
                
                //speed 3 ex=0 ey=1 ez=0
                X_t=X; Y_t=Y+1; Z_t=Z;
                if(Y_t==Ny) Y_t=0;
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[3*Nx*Ny*Nz+tid_t]=f3;
                fOut[getIdx(nnodes, numSpd, tid_t,3)]=f3;
                
                //speed 4 ex=0 ey=-1 ez=0
                X_t=X; Y_t=Y-1; Z_t=Z;
                if(Y_t<0) Y_t=(Ny-1);
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                ///	fOut[4*Nx*Ny*Nz+tid_t]=f4;
                fOut[getIdx(nnodes, numSpd, tid_t,4)]=f4;
                
                
                //speed 5 ex=ey=0 ez=1
                X_t=X; Y_t=Y; Z_t=Z+1;
                //	if(Z_t==Nz) Z_t=0;
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //fOut[5*Nx*Ny*Nz+tid_t]=f5;
                fOut[getIdx(nnodes, numSpd, tid_t,5)]=f5;
                
                //speed 6 ex=ey=0 ez=-1
                X_t=X; Y_t=Y; Z_t=Z-1;
                //	if(Z_t<0) Z_t=(Nz-1);
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[6*Nx*Ny*Nz+tid_t]=f6;
                fOut[getIdx(nnodes, numSpd, tid_t,6)]=f6;
                
                //speed 7 ex=ey=ez=1
                X_t=X+1; Y_t=Y+1; Z_t=Z+1;
                if(X_t==Nx) X_t=0;
                if(Y_t==Ny) Y_t=0;
                //	if(Z_t==Nz) Z_t=0;
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[7*Nx*Ny*Nz+tid_t]=f7;
                fOut[getIdx(nnodes, numSpd, tid_t,7)]=f7;
                
                //speed 8 ex=-1 ey=1 ez=1
                X_t=X-1; Y_t=Y+1; Z_t=Z+1;
                if(X_t<0) X_t=(Nx-1);
                if(Y_t==Ny) Y_t=0;
                //	if(Z_t==Nz) Z_t=0;
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[8*Nx*Ny*Nz+tid_t]=f8;
                fOut[getIdx(nnodes, numSpd, tid_t,8)]=f8;
                
                //speed 9 ex=1 ey=-1 ez=1
                X_t=X+1; Y_t=Y-1; Z_t=Z+1;
                if(X_t==Nx) X_t=0;
                if(Y_t<0) Y_t=(Ny-1);
                //	if(Z_t==Nz) Z_t=0;
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[9*Nx*Ny*Nz+tid_t]=f9;
                fOut[getIdx(nnodes, numSpd, tid_t,9)]=f9;
                
                //speed 10 ex=-1 ey=-1 ez=1
                X_t=X-1; Y_t=Y-1; Z_t=Z+1;
                if(X_t<0) X_t=(Nx-1);
                if(Y_t<0) Y_t=(Ny-1);
                //	if(Z_t==Nz) Z_t=0;
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[10*Nx*Ny*Nz+tid_t]=f10;
                fOut[getIdx(nnodes, numSpd, tid_t,10)]=f10;
                
                //speed 11 ex=1 ey=1 ez=-1
                X_t=X+1; Y_t=Y+1; Z_t=Z-1;
                if(X_t==Nx) X_t=0;
                if(Y_t==Ny) Y_t=0;
                //	if(Z_t<0) Z_t=(Nz-1);
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[11*Nx*Ny*Nz+tid_t]=f11;
                fOut[getIdx(nnodes, numSpd, tid_t,11)]=f11;
                
                //speed 12 ex=-1 ey=1 ez=-1
                X_t=X-1; Y_t=Y+1; Z_t=Z-1;
                if(X_t<0) X_t=(Nx-1);
                if(Y_t==Ny) Y_t=0;
                //	if(Z_t<0) Z_t=(Nz-1);
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[12*Nx*Ny*Nz+tid_t]=f12;
                fOut[getIdx(nnodes, numSpd, tid_t,12)]=f12;
                
                //speed 13 ex=1 ey=-1 ez=-1
                X_t=X+1; Y_t=Y-1; Z_t=Z-1;
                if(X_t==Nx) X_t=0;
                if(Y_t<0) Y_t=(Ny-1);
                //	if(Z_t<0) Z_t=(Nz-1);
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[13*Nx*Ny*Nz+tid_t]=f13;
                fOut[getIdx(nnodes, numSpd, tid_t,13)]=f13;
                
                //speed 14 ex=ey=ez=-1
                X_t=X-1; Y_t=Y-1; Z_t=Z-1;
                if(X_t<0) X_t=(Nx-1);
                if(Y_t<0) Y_t=(Ny-1);
                //	if(Z_t<0) Z_t=(Nz-1);
                tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
                //	fOut[14*Nx*Ny*Nz+tid_t]=f14;
                
                //fOut[tid_t*numSpd+14]=f14;
                fOut[getIdx(nnodes, numSpd, tid_t,14)]=f14;
            }
        }
    }
}

void OpenChannel3D::stream_out_collect(bool isEven,const int z_start,float * buff_out, const int numStreamSpeeds, const int * streamSpeeds){
    int Ny = this->Ny;
    int Nx = this->Nx;
    int numSpd = this->numSpd;
    int nnodes = this->nnodes;
    
    float * fIn_b;
    if(isEven) {
      fIn_b = fOdd;
    }else{
      fIn_b = fEven;
    }
    
    dummyUse(nnodes,numSpd);
    
    #pragma acc parallel loop collapse(3) \
        present(streamSpeeds[0:numStreamSpeeds]) \
        present(fIn_b[0:nnodes*numSpd]) \
        copyout(buff_out[0:Nx*Ny*numStreamSpeeds*HALO])
    for(int z=0;z<HALO;z++){
        for(int y=0;y<Ny;y++){
            for(int x=0;x<Nx;x++){
                for(int spd=0;spd<numStreamSpeeds;spd++){
                    int tid_l = x+y*Nx+z*Nx*Ny; int tid_g = x+y*Nx+(z+z_start)*Nx*Ny;
                    int stream_spd=streamSpeeds[spd];
                    buff_out[tid_l*numStreamSpeeds+spd]=fIn_b[getIdx(nnodes, numSpd, tid_g,stream_spd)];
                }
            }
        }
    }
}

void OpenChannel3D::stream_in_distribute(bool isEven,const int z_start, const float * buff_in, const int numStreamSpeeds, const int * streamSpeeds){
    int Nx = this->Nx;
    int Ny = this->Ny;
    int numSpd = this->numSpd;
    int nnodes = this->nnodes;
    
    float * fIn_b;
    
    if (isEven) {
      fIn_b = fOdd;
    }else{
      fIn_b = fEven;
    }
    
    dummyUse(nnodes,numSpd);
    #pragma acc parallel loop collapse(3) \
        present(streamSpeeds[0:numStreamSpeeds]) \
        present(fIn_b[0:nnodes*numSpd]) \
        copyin(buff_in[0:Nx*Ny*numStreamSpeeds*HALO])
    for(int z=0;z<HALO;z++){
        for(int y=0;y<Ny;y++){
            for(int x=0;x<Nx;x++){
                for(int spd=0;spd<numStreamSpeeds;spd++){
                    int tid_l=x+y*Nx+z*Nx*Ny; int tid_g = x+y*Nx+(z+z_start)*Nx*Ny;
                    int stream_spd=streamSpeeds[spd];
                    fIn_b[getIdx(nnodes, numSpd, tid_g,stream_spd)]=buff_in[tid_l*numStreamSpeeds+spd];
                }
            }
        }
    }
}

void OpenChannel3D::take_lbm_timestep(bool isEven, MPI_Comm comm){
    // collide and stream lower boundary slices
    D3Q15_process_slices(isEven,HALO,HALO+1);
    
    // collect data from lower HALO slice z = 0
    stream_out_collect(isEven,0,ghost_out_m,numMspeeds,Mspeeds);
    // begin communication to ghost_p_in
    MPI_Isend(ghost_out_m,numHALO,MPI_FLOAT,nd_m,tag_d,comm, &rq_out1);
    MPI_Irecv(ghost_in_p,numHALO,MPI_FLOAT,nd_p,tag_d,comm,&rq_in1);
    // collide and stream upper boundary slices
    D3Q15_process_slices(isEven,totalSlices-2*HALO,totalSlices-HALO);
    
    // collect data from upper HALO slice; z = totalSlices-1
       
    stream_out_collect(isEven,totalSlices-HALO,ghost_out_p,numPspeeds,Pspeeds);
    // begin communication to ghost_m_in
    MPI_Isend(ghost_out_p,numHALO,MPI_FLOAT,nd_p,tag_u,comm,&rq_out2);
    MPI_Irecv(ghost_in_m,numHALO,MPI_FLOAT,nd_m,tag_u,comm,&rq_in2);
    // collide and stream interior lattice points
    D3Q15_process_slices(isEven,HALO+1,totalSlices-2*HALO);
    // ensure communication of boundary lattice points is complete
    MPI_Wait(&rq_in1,&stat);
    MPI_Wait(&rq_in2,&stat);
    
    // copy data from i+1 partition into upper boundary slice
      
    stream_in_distribute(isEven,totalSlices-2*HALO,ghost_in_p,numMspeeds,Mspeeds);
    
    
    // copy data from i-1 partition into lower boundary slice
    stream_in_distribute(isEven,HALO,ghost_in_m,numPspeeds,Pspeeds);
}

void OpenChannel3D::initialize_mpi_buffers(){
    // integer for arithmetic into pointer arrays for nodes I need to
    // communicate to neighbors.
    firstNdm = Nx*Ny*HALO;
    lastNdm = Nx*Ny*(HALO+1);
    firstNdp = nnodes-2*(Nx*Ny*HALO);
    lastNdp = nnodes-(Nx*Ny*HALO);
    
    
    
    numHALO = (Nx*Ny*numPspeeds*HALO);
    
    ghost_in_m = new float[Nx*Ny*numMspeeds*HALO];
    ghost_out_m = new float[Nx*Ny*numMspeeds*HALO];
    ghost_in_p = new float[Nx*Ny*numPspeeds*HALO];
    ghost_out_p = new float[Nx*Ny*numPspeeds*HALO];
    offset = firstSlice*Nx*Ny*sizeof(float);
    numEntries = numMySlices*Nx*Ny;
    
    rho_l = new float[numEntries];
    ux_l = new float[numEntries];
    uy_l = new float[numEntries];
    uz_l = new float[numEntries];
}

void OpenChannel3D::initialize_local_partition_variables(){
    numMySlices = Nz/size;
    if(rank<(Nz%size))
    numMySlices+=1;
    
    firstSlice=(Nz/size)*rank;
    if((Nz%size)<rank){
        firstSlice+=(Nz%size);
    }else{
        firstSlice+=rank;
    }
    lastSlice = firstSlice+numMySlices-1;
    totalSlices=numMySlices+2*HALO;
    nnodes = totalSlices*Nx*Ny;
    
    fEven = new float[nnodes*numSpd];
    fOdd = new float[nnodes*numSpd];
    snl = new int[nnodes];
    inl = new int[nnodes];
    onl = new int[nnodes];
    u_bc = new float[nnodes];
    
    // cout << "rank: " << rank << "In initialize_local_partition_variables:" << endl << " fEven = " << fEven << endl;
    // cout << "rank: " << rank << " fOdd = " << fOdd << endl;
    
    nd_m = rank-1;
    if(nd_m<0)
    nd_m=(size-1);
    
    nd_p = rank+1;
    if(nd_p==size)
    nd_p=0;
    
    // initialize the variables
    int tid;
    for(int z=0; z<totalSlices;z++){
        for(int y=0;y<Ny;y++){
            for(int x=0;x<Nx;x++){
                tid=x+y*Nx+z*Nx*Ny;
                for(int spd=0;spd<numSpd;spd++){
                    fEven[getIdx(nnodes, numSpd, tid,spd)]=rho_lbm*w[spd];
                   //fOdd[getIdx(nnodes, numSpd, tid,spd)]=rho_lbm*w[spd];
                }
            }
        }
    }
    
    // initialize all node lists to zero...
    int tid_l, x,z, y;
    for(int nd=0;nd<nnodes;nd++){
        inl[nd]=0;
        onl[nd]=0;
        snl[nd]=0;
    }
    // initialize the solid node list
    y=0;
    for(int z=0;z<totalSlices;z++){
        for(int x=0;x<Nx;x++){
            tid_l=x+y*Nx+z*Nx*Ny;
            snl[tid_l]=1;
        }
    }
    y=(Ny-1);
    for(int z=0;z<totalSlices;z++){
        for(int x=0;x<Nx;x++){
            tid_l=x+y*Nx+z*Nx*Ny;
            snl[tid_l]=1;
        }
    }
    
    // near and far wall
    for(int z=0;z<totalSlices;z++){
        for(int y=0;y<Ny;y++){
            x = 0;
            tid_l = x+y*Nx+z*Nx*Ny;
            snl[tid_l]=1;
            x = (Nx-1);
            tid_l = x+y*Nx+z*Nx*Ny;
            snl[tid_l]=1;
        }
    }
    
    // Due to how the domain is partitioned, the inlet nodes
    // are all assigned to rank 0 process
    if(rank==0){
        z=1; // to account for the HALO nodes on rank 0, this is z=1, not z=0
        for(int y=1;y<(Ny-1);y++){//<-- skip top and bottom
            for(int x=1;x<(Nx-1);x++){
                tid_l = x+y*Nx+z*Nx*Ny;
                inl[tid_l]=1;
            }
        }
    }
    // again, due to the partitioning strategy, all of the outlet nodes
    // are on the (size-1) partition
    if(rank==(size-1)){
        z=totalSlices-1; // again, to account for the HALO at the outlet, this is z = totalSlices-1, not z=totalSlices...
        for(int y=1;y<(Ny-1);y++){
            for(int x=1;x<(Nx-1);x++){
                tid_l=x+y*Nx+z*Nx*Ny;
                inl[tid_l]=1;
            }
        }
    }
    
    // rank 0 partition has the outlet slice on the left halo (due to logical periodicity)
    if(rank==0){
        z=0; // inlet halo is he outlet slice, so z=0 on the rank 0 process is the outlet slice.
        for(int y=1;y<(Ny-1);y++){
            for(int x=1;x<(Nx-1);x++){
                tid_l=x+y*Nx+z*Nx*Ny;
                onl[tid_l]=1;
            }
        }
    }
    // rank (size-1) has the inlet slice on the right halo. (due to logical periodicity)
    if(rank==(size-1)){
        z=totalSlices-2; // because of the HALO, the outlet nodes are at totalSlices-2, not totalSlices-1
        for(int y=1;y<(Ny-1);y++){
            for(int x=1;x<(Nx-1);x++){
                tid_l=x+y*Nx+z*Nx*Ny;
                onl[tid_l]=1;
            }
        }
    }
    
    // initialize u_bc
    float b = ((float)Ny-1.)/2.;
    float h;
    for(int z=0;z<totalSlices;z++){
        for(int y=0;y<Ny;y++){
            for(int x=0;x<Nx;x++){
                tid=x+y*Nx+z*Nx*Ny;
                if((inl[tid]==1)|(onl[tid]==1)){
                    h=((float)y-b)/b;
                    u_bc[tid]=umax_lbm*(1.-(h*h));
                }else{
                    u_bc[tid]=0.;
                }
            }
        }
    }
    
    // initialize outlet density
    rho_out = rho_lbm;
}

void OpenChannel3D::read_input_file(const string input_file){
    ifstream input_params(input_file.c_str(),ios::in);
    if(!input_params.is_open())
    throw std::runtime_error("Could not open params file!");
    input_params >> LatticeType;
    input_params >> Num_ts;
    input_params >> ts_rep_freq;
    input_params >> plot_freq;
    /*input_params >> obst_type;
    input_params >> obst_param1;
    input_params >> obst_param2;
    input_params >> obst_param3;
    input_params >> obst_param4;*/
    input_params >> rho_lbm;
    input_params >> umax_lbm;
    input_params >> omega;
    input_params >> Nx;
    input_params >> Ny;
    input_params >> Nz;
    input_params.close();
}

void OpenChannel3D::initialize_lattice_data(){
    switch(LatticeType){
        case(1):
        //numSpd=15;
        ex = ex15;
        ey = ey15;
        ez = ez15;
        w = w15;
        numPspeeds = numPspeedsD3Q15;
        numMspeeds = numMspeedsD3Q15;
        Mspeeds = MspeedsD3Q15;
        Pspeeds = PspeedsD3Q15; 
    }
}
