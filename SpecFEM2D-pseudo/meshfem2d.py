import numpy as np

def read_material_table():
    #shared/read_material_table.f90
    global icodemat,rho_s_read,cp,cs,comp_g,QKappa,Qmu,phi_read
    mat_table=[1,1,3.4,6,4,0,0,9999,9999,0,0,0,0,0,0]
    icodemat=np.array([mat_table[1]]) #(nbmodels,)
    rho_s_read=np.array([mat_table[2]])
    cp=np.array([mat_table[3]])
    cs=np.array([mat_table[4]])
    comp_g=np.array([mat_table[5]])
    QKappa=np.array([mat_table[7]])
    Qmu=np.array([mat_table[8]])
    phi_read=np.array([0]) #elasti

def read_interfaces_file():
    #shared/read_interfaces_file.f90
    global nxread,nzread
    zmax=300
    xinterface_coords=np.array([[0,0],[xmax,xmax]]) #(max_npoints_interface,number_of_interfaces)
    zinterface_coords=np.array([[0,zmax],[0,zmax]])
    npoints_of_interfaces=np.array([2,2]) #(number_of_interfaces,)
    nz_layer=np.array([60]) #(number_of_layers,), number_of_layers=number_of_interfaces-1
    nxread=nx_param
    nzread=sum(nz_layer)

def read_regions():
    #shared/read_regions.f90
    global num_material
    reg_table=[1,nxread,1,nzread,1]
    imaterial_number=reg_table[-1]
    vpregion=cp[imaterial_number-1]
    vsregion=cs[imaterial_number-1]
    poisson=0.5*(vpregion**2-2*vsregion**2)/(vpregion**2-vsregion**2)
    num_material=imaterial_number #actually (nelmnts,)

def read_parameter_file():
    #shared/read_parameter_file.F90
    global nelmnts
    read_material_table()
    read_interfaces_file()
    nelmnts=nxread*nzread
    read_regions()

    #shared/read_parameter_file.F90/read_parameter_file_derive_flags()
    if pml_boundary_conditions:
        any_abs=True

#meshfem2D/meshfem2D.F90
#shared/parallel.F90/init_mpi()
    #shared/read_parameter_file.F90/read_parameter_file_init()
    #shared/read_parameter_file.F90/open_parameter_file()
        #shared/param_reader.c/param_open()
    #shared/read_parameter_file.F90/read_parameter_file_only()
        #shared/read_parameter_file.F90/check_parameters()

#shared/parallel.F90/world_size()

#shared/parallel.F90/world_rank()
myrank=0

read_parameter_file()

elmnts=np.empty(shape=ngnod*nelmnts)
num_elmnt=0
for j in range(nzread):
    for i in range(nxread):
        elmnts[num_elmnt*ngnod]=j*(nxread+1)+i
        elmnts[num_elmnt*ngnod+1]=j*(nxread+1)+i+1
        elmnts[num_elmnt*ngnod+2]=(j+1)*(nxread+1)+i+1
        elmnts[num_elmnt*ngnod+3]=(j+1)*(nxread+1)+i
        num_elmnt+=1
