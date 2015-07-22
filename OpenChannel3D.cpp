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
    float rho_lbm = this->rho_lbm;
    
    dummyUse(nnodes);
   
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
                float cu,rho,ux,uy,uz,w;
                int X_t,Y_t,Z_t,tid_t,tid;
                
                tid=X+Y*Nx+Z*Nx*Ny;
                
                //load the data into registers
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
                
                	// set macroscopic boundary conditions
		if(inl[tid]==1){
		  ux=0;uy=0; uz=u_bc[tid];
		  //set rho based on uz
		  rho = (1./(1.-uz))*(2.0*(f6+f11+f12+f13+f14)+(f0+f1+f2+f3+f4));
		}
		if(onl[tid]==1){
		  ux=0.; uy=0.; rho=rho_lbm;
		  uz = -1.+((2.*(f5+f7+f8+f9+f10)+(f0+f1+f2+f3+f4)))/rho;
		}
		if(snl[tid]==1){
		  ux=0.; uy=0.; uz=0.;
		}
	
	
		 //everyone compute equilibrium
		float fe0,fe1,fe2,fe3,fe4,fe5,fe6,fe7,fe8,fe9,fe10,fe11,fe12,fe13,fe14;
		//speed 0 ex=ey=ez=0 w=2./9.
	
		fe0=rho*(2./9.)*(1.-1.5*(ux*ux+uy*uy+uz*uz));
	

		//speed 1 ex=1 ey=ez=0 w=1./9.
		cu=3.*(1.*ux);
		fe1=rho*(1./9.)*(1.+cu+0.5*(cu*cu)-
				1.5*(ux*ux+uy*uy+uz*uz));
	      

		//speed 2 ex=-1 ey=ez=0 w=1./9.
		cu=3.*((-1.)*ux);
		fe2=rho*(1./9.)*(1.+cu+0.5*(cu*cu)-
				1.5*(ux*ux+uy*uy+uz*uz));
	      

		//speed 3 ex=0 ey=1 ez=0 w=1./9.
		cu=3.*(1.*uy);
		fe3=rho*(1./9.)*(1.+cu+0.5*(cu*cu)-
				1.5*(ux*ux+uy*uy+uz*uz));
	      

		//speed 4 ex=0 ey=-1 ez=0 w=1./9.
		cu=3.*(-1.*uy);
		fe4=rho*(1./9.)*(1.+cu+0.5*(cu*cu)-
				1.5*(ux*ux+uy*uy+uz*uz));
	

		//speed 5 ex=ey=0 ez=1 w=1./9.
		cu=3.*(1.*uz);
		fe5=rho*(1./9.)*(1.+cu+0.5*(cu*cu)-
				1.5*(ux*ux+uy*uy+uz*uz));
	

		//speed 6 ex=ey=0 ez=-1 w=1./9.
		cu=3.*(-1.*uz);
		fe6=rho*(1./9.)*(1.+cu+0.5*(cu*cu)-
				1.5*(ux*ux+uy*uy+uz*uz));
	
		//speed 7 ex=ey=ez=1 w=1./72.
		cu=3.*(ux+uy+uz);
		fe7=rho*(1./72.)*(1.+cu+0.5*(cu*cu)-
				  1.5*(ux*ux+uy*uy+uz*uz));
	      

		//speed 8 ex=-1 ey=ez=1 w=1./72.
		cu=3.*(-ux+uy+uz);
		fe8=rho*(1./72.)*(1.+cu+0.5*(cu*cu)-
				  1.5*(ux*ux+uy*uy+uz*uz));
	

		//speed 9 ex=1 ey=-1 ez=1 w=1./72.
		cu=3.*(ux-uy+uz);
		fe9=rho*(1./72.)*(1.+cu+0.5*(cu*cu)-
				  1.5*(ux*ux+uy*uy+uz*uz));
	

		//speed 10 ex=-1 ey=-1 ez=1 w=1/72
		cu=3.*(-ux-uy+uz);
		fe10=rho*(1./72.)*(1.+cu+0.5*(cu*cu)-
				  1.5*(ux*ux+uy*uy+uz*uz));
	      

		//speed 11 ex=1 ey=1 ez=-1 w=1/72
		cu=3.*(ux+uy-uz);
		fe11=rho*(1./72.)*(1.+cu+0.5*(cu*cu)-
				  1.5*(ux*ux+uy*uy+uz*uz));
	

		//speed 12 ex=-1 ey=1 ez=-1 w=1/72
		cu=3.*(-ux+uy-uz);
		fe12=rho*(1./72.)*(1.+cu+0.5*(cu*cu)-
				  1.5*(ux*ux+uy*uy+uz*uz));
	      

		//speed 13 ex=1 ey=ez=-1 w=1/72
		cu=3.*(ux-uy-uz);
		fe13=rho*(1./72.)*(1.+cu+0.5*(cu*cu)-
				  1.5*(ux*ux+uy*uy+uz*uz));
	      

		//speed 14 ex=ey=ez=-1 w=1/72
		cu=3.*(-ux-uy-uz);
		fe14=rho*(1./72.)*(1.+cu+0.5*(cu*cu)-
				  1.5*(ux*ux+uy*uy+uz*uz));
	
		// if on inlet or outlet, compute and bounce-back non-equilibrium part of f.
		if((inl[tid]==1)|(onl[tid]==1)){

		  float ft1,ft2,ft3,ft4,ft5,ft6,ft7,ft8,ft9,ft10,ft11,ft12,ft13,ft14;
		  if(inl[tid]==1){
		    //adjust fIn for the unknown velocities: 5,7,8,9,10
		    //bounce-back of non-equilibrium parts
		    //f5, bb_spd=f6
		    f5=fe5+(f6-fe6); //fIn[5*nnodes+tid]=f5;
		    //f7, bb_spd=f14
		    f7=fe7+(f14-fe14); //fIn[7*nnodes+tid]=f7;
		    //f8, bb_spd=f13
		    f8=fe8+(f13-fe13); //fIn[8*nnodes+tid]=f8;
		    //f9, bb_spd=f12
		    f9=fe9+(f12-fe12); //fIn[9*nnodes+tid]=f9;
		    //f10, bb_spd=f11
		    f10=fe10+(f11-fe11); //fIn[10*nnodes+tid]=f10;
		  }else{
		    f6=fe6+(f5-fe5); 
		    f11=fe11+(f10-fe10); 
		    f12=fe12+(f9-fe9); 
		    f13=fe13+(f8-fe8); 
		    f14=fe14+(f7-fe7); 

		  }
		  //ft0=f0-fe0;
		  ft1=f1-fe1; 
		  ft2=f2-fe2;
		  ft3=f3-fe3;
		  ft4=f4-fe4;
		  ft5=f5-fe5;
		  ft6=f6-fe6;
		  ft7=f7-fe7;
		  ft8=f8-fe8;
		  ft9=f9-fe9;
		  ft10=f10-fe10;
		  ft11=f11-fe11;
		  ft12=f12-fe12;
		  ft13=f13-fe13;
		  ft14=f14-fe14;

		  //now, multiply by f# = ((ft#)*Q_flat)*Q_flat'
		  f0= - ft1/3. - ft2/3. - ft3/3. - ft4/3. - ft5/3. - ft6/3. - ft7 - ft8 - ft9 - ft10 - ft11 - ft12 - ft13 - ft14; 
		  f1=(2.*ft1)/3. + (2.*ft2)/3. - ft3/3. - ft4/3. - ft5/3. - ft6/3.; 
		  f2=(2.*ft1)/3. + (2.*ft2)/3. - ft3/3. - ft4/3. - ft5/3. - ft6/3.; 
		  f3=(2.*ft3)/3. - ft2/3. - ft1/3. + (2.*ft4)/3. - ft5/3. - ft6/3.; 
		  f4=(2.*ft3)/3. - ft2/3. - ft1/3. + (2.*ft4)/3. - ft5/3. - ft6/3.; 
		  f5=(2.*ft5)/3. - ft2/3. - ft3/3. - ft4/3. - ft1/3. + (2.*ft6)/3.; 
		  f6=(2.*ft5)/3. - ft2/3. - ft3/3. - ft4/3. - ft1/3. + (2.*ft6)/3.; 
		  f7=(2.*ft1)/3. + (2.*ft2)/3. + (2.*ft3)/3. + (2.*ft4)/3. + (2.*ft5)/3. + (2.*ft6)/3. + 8.*ft7 + 8.*ft14;
		  f8= (2.*ft1)/3. + (2.*ft2)/3. + (2.*ft3)/3. + (2.*ft4)/3. + (2.*ft5)/3. + (2.*ft6)/3. + 8.*ft8 + 8.*ft13;
		  f9= (2.*ft1)/3. + (2.*ft2)/3. + (2.*ft3)/3. + (2.*ft4)/3. + (2.*ft5)/3. + (2.*ft6)/3. + 8.*ft9 + 8.*ft12;
		  f10= (2.*ft1)/3. + (2.*ft2)/3. + (2.*ft3)/3. + (2.*ft4)/3. + (2.*ft5)/3. + (2.*ft6)/3. + 8.*ft10 + 8.*ft11;
		  f11= (2.*ft1)/3. + (2.*ft2)/3. + (2.*ft3)/3. + (2.*ft4)/3. + (2.*ft5)/3. + (2.*ft6)/3. + 8.*ft10 + 8.*ft11;
		  f12= (2.*ft1)/3. + (2.*ft2)/3. + (2.*ft3)/3. + (2.*ft4)/3. + (2.*ft5)/3. + (2.*ft6)/3. + 8.*ft9 + 8.*ft12;
		  f13= (2.*ft1)/3. + (2.*ft2)/3. + (2.*ft3)/3. + (2.*ft4)/3. + (2.*ft5)/3. + (2.*ft6)/3. + 8.*ft8 + 8.*ft13;
		  f14= (2.*ft1)/3. + (2.*ft2)/3. + (2.*ft3)/3. + (2.*ft4)/3. + (2.*ft5)/3. + (2.*ft6)/3. + 8.*ft7 + 8.*ft14;

		  //update fIn for all velocities based on strain tensor
		  //f0, still equals 0..
		  cu = 9./2.; w = 1./9.;

		  //fIn[..] = fe#+f#
		  f0=fe0+f0;

		  f1=fe1+f1*(cu)*w;
		  f2=fe2+f2*(cu)*w;
		  f3=fe3+f3*cu*w;
		  f4=fe4+f4*cu*w;
		  f5=fe5+f5*cu*w;
		  f6=fe6+f6*cu*w;
		  w = 1./72.;
		  f7=fe7+f7*cu*w;
		  f8=fe8+f8*cu*w;
		  f9=fe9+f9*cu*w;
		  f10=fe10+f10*cu*w;
		  f11=fe11+f11*cu*w;
		  f12=fe12+f12*cu*w;
		  f13=fe13+f13*cu*w;
		  f14=fe14+f14*cu*w;
		}

		//everyone relax...
		f0=f0-omega*(f0-fe0);
		f1=f1-omega*(f1-fe1);
		f2=f2-omega*(f2-fe2);
		f3=f3-omega*(f3-fe3);
		f4=f4-omega*(f4-fe4);
		f5=f5-omega*(f5-fe5);
		f6=f6-omega*(f6-fe6);
		f7=f7-omega*(f7-fe7);
		f8=f8-omega*(f8-fe8);
		f9=f9-omega*(f9-fe9);
		f10=f10-omega*(f10-fe10);
		f11=f11-omega*(f11-fe11);
		f12=f12-omega*(f12-fe12);
		f13=f13-omega*(f13-fe13);
		f14=f14-omega*(f14-fe14);
		// stream data
                
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
    //for(int z=0;z<totalSlices;z++){
    //    for(int y=0;y<Ny;y++){
    //        x = 0;
    //        tid_l = x+y*Nx+z*Nx*Ny;
    //        snl[tid_l]=1;
    //        x = (Nx-1);
    //        tid_l = x+y*Nx+z*Nx*Ny;
    //        snl[tid_l]=1;
    //    }
    //}
    
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
    //float b = ((float)Ny-1.)/2.;
    //float h;
    for(int z=0;z<totalSlices;z++){
        for(int y=0;y<Ny;y++){
            for(int x=0;x<Nx;x++){
                tid=x+y*Nx+z*Nx*Ny;
                if((inl[tid]==1)|(onl[tid]==1)){
      //              h=((float)y-b)/b;
      //              u_bc[tid]=umax_lbm*(1.-(h*h));
                  u_bc[tid]=umax_lbm;
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
    input_params >> Warmup_ts;
    input_params >> plot_freq;
    input_params >> Cs;
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
