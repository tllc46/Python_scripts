from obspy import read_inventory
from obspy.core.inventory import PolesZerosResponseStage,CoefficientsTypeResponseStage,ResponseStage
import numpy as np

digital_coeff_norm_tol=0.02 #digital coefficient normalization tolerance
use_estimated_delay=False #추정 지연 사용
use_reported_sensitivity=False #기록된 감도 사용
frequency=10
output_unit="velocity"

def analog_pz_response(frequency):
    #전달 함수 유형 확인
    if stage.pz_transfer_function_type[:7]=="LAPLACE":
        frequency=2*np.pi*frequency
    s=1j*frequency

    #분모 계산
    denominator=1
    for pole in stage.poles:
        denominator*=s-pole

    #분자 계산
    numerator=1
    for zero in stage.zeros:
        numerator*=s-zero

    return numerator/denominator

def digital_coeff_response(frequency):
    num_coeff=len(stage.numerator)
    coefficient=np.array(object=stage.numerator)
    half=(num_coeff-1)/2
    degree=half-np.arange(stop=num_coeff)
    z=np.exp(2*np.pi*1j*frequency*degree/stage.decimation_input_sample_rate)

    return sum(coefficient*z)

#main
inv=read_inventory(path_or_file_object="RESP.NS.N102..HHE",format="RESP",level="response")
response=inv[0][0][0].response
f0=response.instrument_sensitivity.frequency #감도를 결정하는 진동수
response_value=1
calc_sensitivity=1 #계산된 감도
delay=0 #총 지연

for stage in response.response_stages:
    if type(stage)==CoefficientsTypeResponseStage: #digital 계수
        #0Hz에서 정규화 되었는지 확인
        coefficient=np.array(object=stage.numerator)
        summation=sum(coefficient)
        if abs(1-summation)>digital_coeff_norm_tol: #0Hz에서 정규화 되어 있지 않다
            stage.numerator=list(coefficient/summation) #재정규화

        #대칭 유형 결정
        num_coeff=len(stage.numerator) #분자 개수
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
    
    if type(stage)!=ResponseStage: #gain만 있는 단계 제외
        #gain 조정과 재정규화
        if stage.stage_gain_frequency!=f0:
            if type(stage)==PolesZerosResponseStage: #analog 극점과 영점
                R=abs(analog_pz_response(frequency=stage.stage_gain_frequency))
                R0=abs(analog_pz_response(frequency=f0))
                stage.normalization_frequency=f0
            elif type(stage)==CoefficientsTypeResponseStage: #digital 계수
                R=abs(digital_coeff_response(frequency=stage.stage_gain_frequency))
                R0=abs(digital_coeff_response(frequency=f0))
            stage.normalization_factor=1/R0 #재정규화
            stage.stage_gain=stage.stage_gain*R0/R #gain도 조정
            stage.stage_gain_frequency=f0

        #재정규화만
        else:
            if type(stage)==PolesZerosResponseStage and stage.normalization_frequency!=f0: #analog 극점과 영점
                R0=abs(analog_pz_response(frequency=f0))
                stage.normalization_factor=1/R0 #A0 정규화 인자 변경
                stage.normalization_frequency=f0
            elif type(stage)==CoefficientsTypeResponseStage: #digital 계수
                stage.normalization_factor=1

    #기기 응답 값 계산
    if type(stage)==PolesZerosResponseStage: #analog 극점과 영점
        response_value*=analog_pz_response(frequency=frequency)*stage.normalization_factor
    elif type(stage)==CoefficientsTypeResponseStage: #digial 계수
        response_value*=digital_coeff_response(frequency=frequency)*stage.normalization_factor
        #지연
        if stage.symmetry=="asymmetry": #비대칭 계수
            if use_estimated_delay:
                delay+=stage.decimation_delay
            else:
                calc_delay=((num_coeff-1)/2)/stage.decimation_input_sample_rate
                delay+=stage.decimation_correction-calc_delay

    calc_sensitivity*=stage.stage_gain

#감도
if use_reported_sensitivity:
    response_value*=response.instrument_sensitivity.value*np.exp(2*np.pi*1j*frequency*delay)
else:
    response_value*=calc_sensitivity*np.exp(2*np.pi*1j*frequency*delay)

#출력 단위 변환
if output_unit=="displacement":
    response_value*=2*np.pi*1j*frequency
elif output_unit=="acceleration":
    response_value/=2*np.pi*1j*frequency

print(abs(response_value))
print(np.angle(z=response_value,deg=True))
print(response_value)
