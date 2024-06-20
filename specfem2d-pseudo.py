import numpy as np

def compute_coef_convolution(bb):
    #specfem2D/pml_compute.f90
    temp=np.exp(-0.5*bb*dt)
    coef0=temp**2
    coef1=np.empty(shape=(ngllx,ngllz))
    coef2=np.empty(shape=(ngllx,ngllz))
    idx1=1e-5<abs(bb)
    coef1[idx1]=(1-temp[idx1])/bb[idx1]
    coef2[idx1]=coef1[idx1]*temp[idx1]
    idx2=abs(bb)<=1e-5
    coef1[idx2]=0.5*dt
    coef2[idx2]=coef1[idx2]
    return coef0,coef1,coef2

def lik_parameter_computation(kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z,index_ik):
    #specfem2D/pml_compute.f90
    if index_ik==13:
        cpml_x_only_temp=cpml_x_only
        cpml_z_only_temp=cpml_z_only
    elif index_ik==31:
        cpml_x_only_temp=cpml_z_only
        cpml_z_only_temp=cpml_x_only

    if cpml_region_local==cpml_xz:
        a_0=kappa_x/kappa_z
        a_1=-a_0*(alpha_x-alpha_z)*(alpha_x-beta_x)/(alpha_x-beta_z)
        a_2=-a_0*(beta_z-alpha_z)*(beta_z-beta_x)/(beta_z-alpha_x)
    elif cpml_region_local==cpml_x_only_temp:
        a_0=kappa_x
        a_1=-a_0*(alpha_x-beta_x)
        a_2=0
    elif cpml_region_local==cpml_z_only_temp:
        a_0=1/kappa_z
        a_1=0
        a_2=-a_0*(beta_z-alpha_z)

    coef0_1,coef1_1,coef2_1=compute_coef_convolution(alpha_x)
    coef0_2,coef1_2,coef2_2=compute_coef_convolution(beta_z)

    return a_0,a_1,a_2,coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2

def l_parameter_computation():
    #specfem2D/pml_compute.f90
    if cpml_region_local==cpml_xz:
        a_0=kappa_x*kappa_z
        a_1=a_0*((beta_x-alpha_x)+(beta_z-alpha_z))
        a_2=a_0*((beta_x-alpha_x)*(beta_z-alpha_z)-alpha_x*(beta_x-alpha_x)-alpha_z*(beta_z-alpha_z))
        a_3=a_0*alpha_x**2*(alpha_x-beta_x)*(alpha_x-beta_z)/(alpha_z-alpha_x)
        a_4=a_0*alpha_z**2*(alpha_z-beta_x)*(alpha_z-beta_z)/(alpha_x-alpha_z)
    elif cpml_region_local==cpml_x_only:
        a_0=kappa_x
        a_1=a_0*(beta_x-alpha_x)
        a_2=-a_0*alpha_x*(beta_x-alpha_x)
        a_3=a_0*alpha_x**2*(beta_x-alpha_x)
        a_4=0
    elif cpml_region_local==cpml_z_only:
        a_0=kappa_z
        a_1=a_0*(beta_z-alpha_z)
        a_2=-a_0*alpha_z*(beta_z-alpha_z)
        a_3=0
        a_4=a_0*alpha_z**2*(beta_z-alpha_z)

    coef0_1,coef1_1,coef2_1=compute_coef_convolution(alpha_x)
    coef0_2,coef1_2,coef2_2=compute_coef_convolution(alpha_z)

    return a_0,a_1,a_2,a_3,a_4,coef0_1,coef1_1,coef2_1,coef0_1,coef1_2,coef2_2

def pml_boundary_elastic():
    #specfem2D/pml_compute.f90
    global displ_elastic_old,displ_elastic,veloc_elastic,accel_elastic

    if not anyabs:
        return

    if not pml_nglob_abs_elastic:
        return

    idx=np.arange(stop=pml_nglob_abs_elastic) #... or maybe idx=range(pml_nglob_abs_elastic)
    iglob=pml_abs_points_elastic(idx)

    displ_elastic_old[:,iglob]=np.zeros(shape=(ndim,pml_nglob_abs_elastic))
    displ_elastic[:,iglob]=np.zeros(shape=(ndim,pml_nglob_abs_elastic))
    veloc_elastic[:,iglob]=np.zeros(shape=(ndim,pml_nglob_abs_elastic))
    accel_elastic[:,iglob]=np.zeros(shape=(ndim,pml_nglob_abs_elastic))

def pml_compute_memory_variables_elastic():
    #specfem2D/pml_compute_memory_variables.f90
    global dux_dxl,dux_dzl,duz_dxl,duz_dzl
    global pml_dux_dxl,pml_dux_dzl,pml_duz_dxl,pml_duz_dzl
    global kappa_x,kappa_z,d_x,d_z,alpha_x,alpha_z,beta_x,beta_z
    global rmemory_dux_dx,rmemory_duz_dx,rmemory_dux_dz,rmemory_duz_dz

    if nspec_pml:
        return

    pml_dux_dxl=dux_dxl
    pml_dux_dzl=dux_dzl
    pml_duz_dzl=duz_dzl
    pml_duz_dxl=duz_dxl

    dux_dxi_old=np.sum(a=displ_elastic_old[0,iglob[None]]*hprime_xx[:,:,None]],axis=1)
    duz_dxi_old=np.sum(a=displ_elastic_old[1,iglob[None]]*hprime_xx[:,:,None]],axis=1)
    dux_dgamma_old=np.sum(a=displ_elastic_old[0,iglob[:,None]]*hprime_zz[None],axis=2)
    duz_dgamma_old=np.sum(a=displ_elastic_old[1,iglob[:,None]]*hprime_zz[None],axis=2)

    pml_dux_dxl_old=dux_dxi_old*xix[:,:,ispec]+dux_dgamma_old*gammax[:,:,ispec]
    pml_dux_dzl_old=dux_dxi_old*xiz[:,:,ispec]+dux_dgamma_old*gammaz[:,:,ispec]
    pml_duz_dxl_old=duz_dxi_old*xix[:,:,ispec]+duz_dgamma_old*gammax[:,:,ispec]
    pml_duz_dzl_old=duz_dxi_old*xiz[:,:,ispec]+duz_dgamma_old*gammaz[:,:,ispec]

    ispec_pml=spec_to_pml(ispec)
    cpml_region_local=region_cpml(ispec)

    kappa_x=K_x_store[:,:,ispec_pml]
    kappa_z=K_z_store[:,:,ispec_pml]

    d_x=d_x_store[:,:,ispec_pml]
    d_z=d_z_store[:,:,ispec_pml]

    alpha_x=alpha_x_store[:,:,ispec_pml]
    alpha_z=alpha_z_store[:,:,ispec_pml]

    beta_x=alpha_x+d_x/kappa_x
    beta_z=alpha_z+d_z/kappa_z

    a5,a6,a7,coef0_zx_1,coef1_zx_1,coef2_zx_1,coef0_zx_2,coef1_zx_2,coef2_zx_2=lik_parameter_computation(kappa_z,beta_z,alpha_z,kappa_x,beta_x,alpha_x,index_ik=31)
    a8,a9,a10,coef0_xz_1,coef1_xz_1,coef2_xz_1,coef0_xz_2,coef1_xz_2,coef2_xz_2=lik_parameter_computation(kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z,index_ik=13)

    rmemory_dux_dx[:,:,ispec_pml,0]=coef0_zx_1*rmemory_dux_dx[:,:,ispec_pml,0]+coef1_zx_1*pml_dux_dxl+coef2_zx_1*pml_dux_dxl_old
    rmemory_duz_dx[:,:,ispec_pml,0]=coef0_zx_1*rmemory_duz_dx[:,:,ispec_pml,0]+coef1_zx_1*pml_duz_dxl+coef2_zx_1*pml_duz_dxl_old
    rmemory_dux_dx[:,:,ispec_pml,1]=coef0_zx_2*rmemory_dux_dx[:,:,ispec_pml,1]+coef1_zx_2*pml_dux_dxl+coef2_zx_2*pml_dux_dxl_old
    rmemory_duz_dx[:,:,ispec_pml,1]=coef0_zx_2*rmemory_duz_dx[:,:,ispec_pml,1]+coef1_zx_2*pml_duz_dxl+coef2_zx_2*pml_duz_dxl_old

    rmemory_dux_dz[:,:,isepc_pml,0]=coef0_xz_1*rmemory_dux_dz[:,:,ispec_pml,0]+coef1_xz_1*pml_dux_dzl+coef2_xz_1*pml_dux_dzl_old
    rmemory_duz_dz[:,:,ispec_pml,0]=coef0_xz_1*rmemory_duz_dz[:,:,ispec_pml,0]+coef1_xz_1*pml_duz_dzl+coef2_xz_1*pml_duz_dzl_old
    rmemory_dux_dz[:,:,ispec_pml,1]=coef0_xz_2*rmemory_dux_dz[:,:,ispec_pml,1]+coef1_xz_2*pml_dux_dzl+coef2_xz_2*pml_dux_dzl_old
    rmemory_dux_dz[:,:,isepc_pml,1]=coef0_xz_2*rmemory_duz_dz[:,:,ispec_pml,1]+coef1_xz_2*pml_duz_dzl+coef2_xz_2*pml_duz_dzl_old

    dux_dxl=a5*pml_dux_dxl+a6*rmemory_dux_dx[:,:,isepc_pml,0]+a7*rmemory_dux_dx[:,:,ispec_pml,1]
    duz_dxl=a5*pml_duz_dxl+a6*rmemory_duz_dx[:,:,ispec_pml,0]+a7*rmemory_duz_dx[:,:,ispec_pml,1]
    dux_dzl=a8*pml_dux_dzl+a9*rmemory_dux_dz[:,:,isepc_pml,0]+a10*rmemory_dux_dz[:,:,ispec_pml,1]
    duz_dzl=a8*pml_duz_dzl+a9*rmemory_duz_dz[:,:,isepc_pml,0]+a10*rmemory_duz_dz[:,:,ispec_pml,1]

def pml_compute_accel_contribution_elastic():
    #specfem2D/pml_compute_accel_contribution.f90
    global accel_elastic_pml
    global rmemory_displ_elastic

    ispec_pml=spce_to_pml(ispec)
    cpml_region_local=region_cpml(ispec)

    fac=rhol*jacobian[:,:,ispec]

    a0,a1,a2,a3,a4,coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2=l_parameter_computation()

    rmemory_displ_elastic[0,0,:,:,ispec_pml]=coef0_1*rmemory_displ_elastic[0,0,:,:,ispec_pml]+coef1_1*dummy_loc[0]+coef2_1*displ_elastic_old[0,iglob]
    rmemory_displ_elastic[0,1,:,:,ispec_pml]=coef0_1*rmemory_displ_elastic[0,1,:,:,ispec_pml]+coef1_1*dummy_loc[1]+coef2_1*displ_elastic_old[1,iglob]

    rmemory_displ_elastic[1,0,:,:,ispec_pml]=coef0_2*rmemory_displ_elastic[1,0,:,:,ispec_pml]+coef1_2*dummy_loc[0]+coef2_2*displ_elastic_old[0,iglob]
    rmemory_displ_elastic[1,1,:,:,ispec_pml]=coef0_2*rmemory_displ_elastic[1,1,:,:,ispec_pml]+coef1_2*dummy_loc[1]+coef2_2*displ_elastic_old[1,iglob]

    accel_elastic_pml[0]=wxgll[:,None]*wzgll[None]*fac*(a1*veloc_elastic[0,iglob]+a2*dummy_loc[0]+a3*rmemory_displ_elastic[0,0,:,:,ispec_pml]+a4*rmemory_displ_elastic[1,0,:,:,ispec_pml])
    accel_elastic_pml[1]=wxgll[:,None]*wzgll[None]*fac*(a1*veloc_elastic[1,iglob]+a2*dummy_loc[1]+a3*rmemory_displ_elastic[0,1,:,:,ispec_pml]+a4*rmemory_displ_elastic[1,1,:,:,ispec_pml])

def compute_forces_viscoelastic():
    #specfem2D/compute_forces_viscoelastic.F90
    global dux_dxl,dux_dzl,duz_dxl,duz_dzl
    global accel_elastic
    global dummy_loc
    global rhol

    for ispec_p in range(num_elements):
        ispec=phase_ispec_inner_elastic(ispec_p,iphase)

        iglob=ibool[:,:,ispec]
        dummy_loc=displ_elastic[:,iglob]

        mxm_4comp_singleA(...)

        dux_dxl=dux_dxi*xix[:,:,ispec]+dux_dgamma*gammax[:,:,ispec]
        dux_dzl=dux_dxi*xiz[:,:,ispec]+dux_dgamma*gammaz[:,:,ispec]

        duz_dxl=duz_dxi*xix[:,:,ispec]+duz_dgamma*gammax[:,:,ispec]
        duz_dzl=duz_dxi*xiz[:,:,ispec]+duz_dgamma*gammaz[:,:,ispec]

        if pml_boundary_conditions and ispec_is_pml(ispec):
            pml_compute_memory_variables_elastic()

        mul_unrelaxed_elastic=mustore[:,:,ispec]
        rhol=rhostore[:,:,ispec]
        cpl=rho_vpstore[:,:,ispec]/rhol

        lambdal_unrelaxed_elastic=rhol*cpl**2-2*mul_unrelaxed_elastic
        lambdaplus2mu_unrelaxed_elastic=lambdal_unrelaxed_elastic+2*mul_unrelaxed_elastic
        lambdaplusmul_unrelaxed_elastic=lambdal_unrelaxed_elastic+mul_unrelaxed_elastic

        simga_xx=lambdaplus2mu_unrelaxed_elastic*dux_dxl+lambdal_unrelaxed_elastic*duz_dzl
        sigma_zz=lambdaplus2mu_unrelaxed_elastic*duz_dzl+lambdal_unrelaxed_elastic*dux_dxl
        sigma_xz=mul_unrelaxed_elastic*(duz_dxl+dux_dzl)
        sigma_zx=sigma_xz

        if pml_boundary_conditions and ispec_is_pml(ispec) and 0<nspec_pml:
            sigma_xx=lambdaplus2mu_unrelaxed_elastic*dux_dxl+lambdal_unrelaxed_elastic*pml_duz_dzl
            sigma_zz=lambdaplus2mu_unrelaxed_elastic*duz_dzl+lambdal_unrelaxed_elastic*pml_dux_dxl
            sigma_zx=mul_unrelaxed_elastic*(pml_duz_dxl+dux_dzl)
            simga_xz=mul_unrelaxed_elastic*(pml_dux_dzl+duz_dxl)

        tempx1=jacobian[:,:,ispec]*(sigma_xx*xix[:,:,ispec]+sigma_zx*xiz[:,:,ispec])
        tempz1=jacobian[:,:,ispec]*(sigma_xz*xix[:,:,ispec]+sigma_zz*xiz[:,:,ispec])

        tempx2=jacobian[:,:,ispec]*(sigma_xx*gammax[:,:,ispec]+sigma_zx*gammaz[:,:,ispec])
        tempz2=jacobian[:,:,ispec]*(sigma_xz*gammax[:,:,ispec]+sigma_zz*gammaz[:,:,ispec])

        if pml_boundary_conditions:
            pml_compute_accel_contribution_elastic()

        if not iglob_is_forced(iglob):
            tempx1l=np.sum(a=tempx1[:,None]*hprimewgll_xx[:,:,None],axis=0)
            tempx2l=np.sum(a=tempx2[:,:,None]*hprimewgll_zz[None],axis=1)
            tempz1l=np.sum(a=tempz1[:,None]*hprimewgll_xx[:,:,None],axis=0)
            tempz2l=np.sum(a=tempz2[:,:,None]*hprimewgll_zz[None],axis=1)

            accel_elastic[0,iglob]-=wzgll[None]*tempx1l+wxgll[:,None]*tempx2l
            accel_elastic[1,iglob]-=wzgll[None]*tempz1l+wxgll[:,None]*tempz2l

        if pml_boundary_conditions and ispec_is_pml(ispec) and not iglob_is_forced(iglob):
            accel_elastic[:,iglob]-=accel_elastic_pml

def update_displ_Newmark():
    #specfem2D/update_displacement_Newmark.F90
    update_displ_elastic_forward()

def update_displ_elastic_forward():
    #specfem2D/update_displacement_Newmark.F90
    update_displacement_newmark_elastic()

def update_displacement_newmark_elastic():
    #specfem2D/update_displacement_Newmark.F90
    global displ_elastic_old,displ_elastic,veloc_elastic,accel_elastic
    if pml_boundary_conditions:
        displ_elastic_old=displ_elastic+deltatsquareover2/2*accel_elastic
    displ_elastic+=deltat*veloc_elastic+deltatsquareover2*accel_elastic
    veloc_elastic+=deltatover2*accel_elastic
    accel_elastic=np.zeros(shape=(ndim,nglob_elastic))

def update_veloc_elastic_Newmark():
    #specfem2D/update_displacement_Newmark.F90
    global veloc_elastic
    veloc_elastic+=deltatover2*accel_elastic

def compute_forces_viscoelastic_main():
    #specfem2D/compute_forces_poroelastic_calling_routine.F90
    global accel_elastic

    for iphase in range(2):
        compute_forces_viscoelastic()

        if iphase==1:
            if pml_boundary_conditions:
                pml_boundary_elastic()

    if pml_boundary_conditions and 0<nglob_interface and save_forward:
        for i in range(nglob_interface):
            print(accel_elastic[:,point_interface[i]],veloc_elastic[:,point_interface[i]],displ_elastic[:,point_interface[i]])

    accel_elastic*=rmass_inverse_elastic

    update_veloc_elastic_Newmark()

def iterate_time():
    #specfem2D/iterate_time.F90
    for it in range(nstep):
        for i_stage in range(nstage_time_scheme):
            update_displ_Newmark()
            compute_forces_viscoelastic_main()
