from obspy import read_inventory
from obspy.core.inventory import PolesZerosResponseStage,CoefficientsTypeResponseStage,ResponseStage
import numpy as np

digital_coeff_norm_tol=0.02 #digital coefficient normalization tolerance
use_estimated_delay=False #use estimated delay 
use_reported_sensitivity=False #use reported sensitivity
frequency=10
output_unit="velocity"

def analog_pz_response(stage,frequency):
    #check transfer function type
    if stage.pz_transfer_function_type[:7]=="LAPLACE":
        frequency=2*np.pi*frequency
    s=1j*frequency

    #calculate denominator
    denominator=1
    for pole in stage.poles:
        denominator*=s-pole

    #calculate numerator
    numerator=1
    for zero in stage.zeros:
        numerator*=s-zero

    return numerator/denominator

def digital_coeff_response(stage,frequency):
    num_coeff=len(stage.numerator)
    coefficient=np.array(stage.numerator)
    half=(num_coeff-1)/2
    degree=np.arange(half,-(half+1),-1)
    z=np.exp(2*np.pi*1j*frequency*degree/stage.decimation_input_sample_rate)

    return sum(coefficient*z)

#main
inv=read_inventory(path_or_file_object="RESP.NS.N102..HHE",format="RESP",level="response")
response=inv[0][0][0].response
f0=response.instrument_sensitivity.frequency #frequency of sensitivity
response_value=1
calc_sensitivity=1 #calculated sensitivity
delay=0 #total delay

for stage in response.response_stages:
    if type(stage)==CoefficientsTypeResponseStage: #digital coefficients
        #check normalization at 0Hz
        coefficient=np.array(stage.numerator)
        summation=sum(coefficient)
        if abs(1-summation)>digital_coeff_norm_tol: #not normlized at 0Hz
            stage.numerator=list(coefficient/summation) #normalization
        
        #determine symmetry type
        num_coeff=len(stage.numerator) #number of numerator
        if num_coeff%2==0:
            half=int(num_coeff/2)
            if stage.numerator[:half]==stage.numerator[half:][::-1]:
                stage.symmetry="even symmetry"
        else:
            half=int((num_coeff-1)/2)
            if stage.numerator[:half]==stage.numerator[half+1:][::-1]:
                stage.symmetry="odd symmetry"
        if not hasattr(stage,"symmetry"):
            stage.symmetry="asymmetry"
    
    if type(stage)!=ResponseStage: #except gain stage
        #gain scaling and normalize again
        if stage.stage_gain_frequency!=f0:
            if type(stage)==PolesZerosResponseStage: #analog poles and zeros
                R=abs(analog_pz_response(stage,stage.stage_gain_frequency))
                R0=abs(analog_pz_response(stage,f0))
                stage.normalization_frequency=f0
            elif type(stage)==CoefficientsTypeResponseStage: #digital coefficients
                R=abs(digital_coeff_response(stage,stage.stage_gain_frequency))
                R0=abs(digital_coeff_response(stage,f0))
            stage.normalization_factor=1/R0 #change normalization
            stage.stage_gain=stage.stage_gain*R0/R #change gain
            stage.stage_gain_frequency=f0

        #only normalize again
        else:
            if type(stage)==PolesZerosResponseStage and stage.normalization_frequency!=f0: #analog poles and zeros
                R0=abs(analog_pz_response(stage,f0))
                stage.normalization_factor=1/R0 #change A0 normalization factor
                stage.normalization_frequency=f0
            elif type(stage)==CoefficientsTypeResponseStage: #digital coefficients
                stage.normalization_factor=1

    #calculate response value
    if type(stage)==PolesZerosResponseStage: #analog poles and zeros
        response_value*=analog_pz_response(stage,frequency)*stage.normalization_factor
    elif type(stage)==CoefficientsTypeResponseStage: #digial coefficients
        response_value*=digital_coeff_response(stage,frequency)*stage.normalization_factor
        #delay
        if stage.symmetry=="asymmetry": #asymmetry type
            if use_estimated_delay:
                delay+=stage.decimation_delay
            else:
                calc_delay=((num_coeff-1)/2)/stage.decimation_input_sample_rate
                delay+=stage.decimation_correction-calc_delay

    calc_sensitivity*=stage.stage_gain

#sensitivity
if use_reported_sensitivity:
    response_value*=response.instrument_sensitivity.value*np.exp(2*np.pi*1j*frequency*delay)
else:
    response_value*=calc_sensitivity*np.exp(2*np.pi*1j*frequency*delay)

#unit conversions
if output_unit=="displacement":
    response_value*=2*np.pi*1j*frequency
elif output_unit=="acceleration":
    response_value/=2*np.pi*1j*frequency

print(abs(response_value))
print(np.angle(response_value,deg=True))
print(response_value)
