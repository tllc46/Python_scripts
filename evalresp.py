from sys import stderr

from obspy import read_inventory
from obspy.core.inventory import PolesZerosResponseStage,CoefficientsTypeResponseStage,ResponseStage
import numpy as np

#private.h
FIR_NORM_TOL=0.02

#private.h/filt_types
FIR_SYM_1=4 #odd number of weights
FIR_SYM_2=5 #even number of weights
FIR_ASYM=6

#private.h/units
DIS=1 #displacement
VEL=2 #velocity
ACC=3 #acceleration
DEFAULT=6 #default

#public.h
REVNUM="5.0.1"

def read_channel_header():
    #input.c/read_channel_header()
    response.calc_delay=0
    response.estim_delay=0
    response.applied_corr=0

def read_channel_data():
    #input.c/read_channel_data()
    for stage in response_stages:
        #alloc_fctns.c/alloc_fir()
        if type(stage)==CoefficientsTypeResponseStage: #blockette 54
            stage.h0=1

        #alloc_fctns.c/alloc_deci()
        stage.sample_int=0

        #input.c/read_deci()
        if stage.decimation_input_sample_rate:
            stage.sample_int=1/stage.decimation_input_sample_rate

def check_symmetry(stage):
    #resp_fncts.c/check_symmetry()
    nc=len(stage.numerator)
    num_sum=sum(stage.numerator)
    if nc and FIR_NORM_TOL<abs(num_sum-1):
        print(f" WARNING: FIR normalized: sum[coef]={num_sum:E}; ",file=stderr)
        print(f" {net.code} {sta.code} {loc} {cha.code}",file=stderr)
        for k in range(nc):
            stage.numerator[k]/=num_sum

    stage.symmetry=FIR_ASYM

    if not nc%2:
        n0=nc//2
        if stage.numerator[:n0]==stage.numerator[n0:][::-1]:
            stage.symmetry=FIR_SYM_2

    else:
        n0=(nc-1)//2
        if stage.numerator[:n0]==stage.numerator[n0+1:][::-1]:
            stage.symmetry=FIR_SYM_1

def check_channel():
    #resp_fncts/check_channel()
    prev_stage=None
    for stage in response_stages:
        deci_flag=False
        if stage.decimation_factor:
            if type(stage)==CoefficientsTypeResponseStage:
                check_symmetry(stage=stage)
                nc=len(stage.numerator)
                response.calc_delay+=(nc-1)/2*stage.sample_int

            response.estim_delay+=stage.decimation_delay
            response.applied_corr+=stage.decimation_correction
            response.sint=stage.sample_int*stage.decimation_factor

            deci_flag=True

        if type(stage)!=ResponseStage: #GAIN
            if not stage.stage_gain:
                print(f"{net.code}.{sta.code}.{loc}.{cha.code}: Missing gain in stage {stage.stage_sequence_number} - using unit gain at 1Hz",file=stderr)
                stage.stage_gain=1
                stage.stage_gain_frequency=1

            if prev_stage and prev_stage.output_units!=stage.input_units:
                print(f"{net.code}.{sta.code}.{loc}.{cha.code}: units mismatch between stages",file=stderr)
                exit()
            prev_stage=stage

        if type(stage)==CoefficientsTypeResponseStage and not deci_flag:
            print(f"{net.code}.{sta.code}.{loc}.{cha.code}: required decimation blockette for IIR or FIR filter missing",file=stderr)
            exit()

def analog_trans(stage,freq):
    #calc_fctns.c/analog_trans()
    omega=1j*freq
    if stage.pz_transfer_function_type=="LAPLACE (RADIANS/SECOND)":
        omega=2*np.pi*1j*freq
    denom=1
    num=1
    h0=stage.normalization_factor

    for ze in stage.zeros:
        num*=omega-ze

    for po in stage.poles:
        denom*=omega-po

    out=h0*num/denom
    return out

def fir_trans(stage,w):
    #calc_fctns.c/fir_sym_trans() and calc_fctns.c/fir_asym_trans()
    a=np.array(object=stage.numerator)
    na=len(a)
    h0=stage.h0
    sint=stage.sample_int #sample interval
    wsint=w*sint

    half=(na-1)/2
    degree=half-np.arange(stop=na)
    z=np.exp(wsint*1j*degree)
    out=h0*sum(a*z)
    return out

def convert_to_units(w):
    #calc_fctns.c/convert_to_units()
    if output_unit==DEFAULT:
        return 1
    elif output_unit not in [DIS,VEL,ACC,DEFAULT]:
        print("convert_to_units: bad output units",file=stderr)
        exit()

    scale_val=1
    inp=response.instrument_sensitivity.input_units
    #first, convert all to M/S
    if inp=="M": #DIS
        if output_unit==DIS:
            return scale_val
        if not w:
            scale_val=-1j/w
        else:
            scale_val=0

    elif inp=="M/S**2": #ACC
        if output_unit==ACC:
            return scale_val
        scale_val=1j*w

    #then convert to final unit
    if output_unit==DIS:
        scale_val*=1j*w

    elif output_unit==ACC:
        if not w:
            scale_val*=-1j/w
        else:
            scale_val=0

    return scale_val

def normalize_response():
    #calc_fctns.c/normalize_response()
    nstages=len(response_stages)

    #test 1
    if nstages==2:
        stage=response_stages[0]
        if not stage.stage_gain:
            if not response.instrument_sensitivity.value:
                print("norm_resp; no stage gain defined, zero sensitivity",file=stderr)
                exit()
            else:
                stage.stage_gain=response.instrument_sensitivity.value
                stage.stage_gain_frequency=response.instrument_sensitivity.frequency

    #test 2
    for stage in response_stages:
        if not stage.stage_gain:
            print("norm_resp; zero stage gain",file=stderr)
            exit()

    if not response.instrument_sensitivity.value:
        for stage in response_stages:
            if stage.stage_gain_frequency:
                response.instrument_sensitivity.frequency=stage.stage_gain_frequency

    calc_sensit=1
    f=response.instrument_sensitivity.frequency
    w=2*np.pi*f #angular frequency omega
    for stage in response_stages:
        if stage.stage_gain:
            if stage.stage_gain_frequency!=f or type(stage)==PolesZerosResponseStage and stage.normalization_frequency!=f:
                reset_gain=True
                if type(stage)==PolesZerosResponseStage: #ANALOG_PZ || LAPLACE_PZ
                    stage.normalization_factor=1
                    df=analog_trans(stage=stage,freq=stage.stage_gain_frequency)
                    if not df:
                        print("norm_resp: Gain frequency of zero found in bandpass analog filter",file=stderr)
                        exit()
                    of=analog_trans(stage=stage,freq=f)
                    if not of:
                        print("norm_resp: Chan. Sens. frequency found with bandpass analog filter",file=stderr)
                        exit()
                elif type(stage)==CoefficientsTypeResponseStage: #FIR_SYM_1 || FIR_SYM_2 || FIR_ASYM
                    df=fir_trans(stage=stage,w=2*np.pi*stage.stage_gain_frequency)
                    of=fir_trans(stage=stage,w=w)
                else:
                    reset_gain=False

                if reset_gain:
                    stage.stage_gain*=abs(of)/abs(df)
                    stage.stage_gain_frequency=f
                    if type(stage)==PolesZerosResponseStage: #ANALOG_PZ || LAPLACE_PZ || IIR_PZ
                        stage.normalization_factor=1/abs(of)
                        stage.normalization_frequency=f
                    elif type(stage)==CoefficientsTypeResponseStage: #FIR_SYM_1 || FIR_SYM_2 || FIR_ASYM
                        stage.h0=1/abs(of)

            calc_sensit*=stage.stage_gain
            if nstages==1:
                response.instrument_sensitivity.value=calc_sensit

    response.calc_sensit=calc_sensit
    #test 3
    if not response.instrument_sensitivity.value:
        percent_diff=abs((response.instrument_sensitivity.value-calc_sensit)/response.instrument_sensitivity.value)
        if 0.05<=percent_diff:
            print(" (norm_resp): computed and reported sensitivities differ by more than 5 percent. ",file=stderr)
            print("\t Execution continuing.",file=stderr)

def calculate_response():
    #calc_fctns.c/calculate_response()
    w=2*np.pi*freq
    val=1
    for stage in response_stages:
        if type(stage)==PolesZerosResponseStage: #ANALOG_PZ || LAPLACE_PZ
            of=analog_trans(stage=stage,freq=freq)

        elif type(stage)==CoefficientsTypeResponseStage: #FIR_SYM_1 || FIR_SYM_2 || FIR_ASYM 
            nc=len(stage.numerator)
            of=fir_trans(stage=stage,w=w)
            estim_delay=stage.decimation_delay
            corr_applied=stage.decimation_correction
            calc_delay=(nc-1)/2*stage.sample_int
            if stage.symmetry==FIR_ASYM:
                if use_estimated_delay:
                    delay=estim_delay
                else:
                    delay=corr_applied-calc_delay
                of*=np.exp(w*delay*1j) #calc_fcnts.c/calc_time_shift()

        val*=of

    if not use_total_sensitivity:
        output=val*response.calc_sensit
    else:
        output=val*response.instrument_sensitivity.value

    output*=convert_to_units(w=w)

    return output

def evalresp_unit_string(output_unit):
    #evaluation.c/evalresp_unit_string()
    if output_unit==DIS:
        return "Displacement (m)"
    elif output_unit==VEL:
        return "Velocity (m/s)"
    elif output_unit==ACC:
        return "Acceleration (m/s**2)"
    elif output_unit==DEF:
        return "Documented response unit"
    else:
        return "Unrecognized response unit"

def evalresp_channel_to_log():
    #output.c/evalresp_channel_to_log()
    print(" --------------------------------------------------",file=stderr)
    print(f"  {filename}",file=stderr)
    print(" --------------------------------------------------",file=stderr)
    print(f"  {net.code} {sta.code} {loc} {cha.code} ",file=stderr)
    print(f"  {cha.start_date} {cha.end_date}",file=stderr)
    print(f"   documented input units: {response.instrument_sensitivity.input_units} - {response.instrument_sensitivity.input_units_description}",file=stderr)
    print(f"   documented output units: {response.instrument_sensitivity.output_units} - {response.instrument_sensitivity.output_units_description}",file=stderr)
    print(f"   requested units: {evalresp_unit_string(output_unit)}",file=stderr)
    print(f"   computed sens={response.calc_sensit:.5E} (reported={response.instrument_sensitivity.value:.5E}) @ {response.instrument_sensitivity.frequency:.5E} Hz",file=stderr)
    print(f"   calc_del={response.calc_delay:.5E}  corr_app={response.applied_corr:.5E}  est_delay={response.estim_delay:.5E}  final_sint={response.sint:.3g}(sec/sample)",file=stderr)
    if use_total_sensitivity:
        print("   (reported sensitivity was used to compute response (-ts option enabled))",file=stderr)

    for stage in response_stages:
        out_str=f"     stage {stage.stage_sequence_number:2}:"
        if type(stage)==PolesZerosResponseStage and stage.pz_transfer_function_type=="LAPLACE (RADIANS/SECOND)": #LAPLACE_PZ
            out_str+=f" LAPLACE     A0={stage.normalization_factor:E} NZeros= {len(stage.zeros):2} NPoles= {len(stage.poles):2}"
        elif type(stage)==PolesZerosResponseStage and stage.pz_transfer_function_type=="LAPLACE (HERTZ)": #ANALOG_PZ
            out_str+=f" ANALOG      A0={stage.normalization_factor:E} NZeros= {len(stage.zeros):2} NPoles= {len(stage.poles):2}"
        elif type(stage)==CoefficientsTypeResponseStage and stage.symmetry==FIR_SYM_1:
            out_str+=f" FIR_SYM_1   H0={stage.h0:E} Ncoeff={len(stage.numerator):3}"
        elif type(stage)==CoefficientsTypeResponseStage and stage.symmetry==FIR_SYM_2:
            out_str+=f" FIR_SYM_2   H0={stage.h0:E} Ncoeff={len(stage.numerator):3}"
        elif type(stage)==CoefficientsTypeResponseStage and stage.symmetry==FIR_ASYM:
            out_str+=f" FIR_ASYM    H0={stage.h0:E} Ncoeff={len(stage.numerator):3}"
        elif type(stage)==ResponseStage: #GAIN
            out_str+=f" GAIN        Sd={stage.stage_gain:E}"
            print(out_str,file=stderr)
            continue

        if stage.decimation_factor:
            out_str+=f" SamInt={stage.sample_int:E}" #DECIMATION
            if stage.decimation_correction<0:
                print(f" WARNING Stage {stage.stage_sequence_number}: Negative correction_applied={stage.decimation_correction:.5E} is likely to be incorrect",file=stderr)
            if stage.decimation_delay<0:
                print(f" WARNING Stage {stage.stage_sequence_number}: Negative estimated_delay={stage.decimation_delay:.5E} is likely to be incorrect",file=stderr)

        out_str+=f" Sd={stage.stage_gain:E}" #GAIN
        print(out_str,file=stderr)

    print("--------------------------------------------------",file=stderr)

#main
filename="RESP.HL.H28..HHZ"
digital_coeff_norm_tol=0.02 #digital coefficient normalization tolerance
use_estimated_delay=False #추정 지연 사용
use_total_sensitivity=False #기록된 감도 사용
freq=10
output_unit=VEL

inv=read_inventory(path_or_file_object=filename,format="RESP",level="response")
net=inv[0]
sta=net[0]
cha=sta[0]
response=cha.response
response_stages=response.response_stages
if not cha.location_code:
    loc="--"
else:
    loc=cha.location_code

print(f"<< EVALRESP RESPONSE OUTPUT V{REVNUM} >>",file=stderr)
read_channel_header()
read_channel_data()
check_channel()
normalize_response()
output=calculate_response()
evalresp_channel_to_log()
