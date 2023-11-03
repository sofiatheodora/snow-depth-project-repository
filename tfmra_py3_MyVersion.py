#author: sofia thirslund
#based on: tfmra_py3 by Alessandro Di Bella provided by Henriette Skourup
import numpy as np


def find_first_peak(wf, th_Pmax):
    """
    Find the first peak above th_Pmax
    """

    for i in range(1, len(wf)):
        # Check on Pmax
        if wf[i] <= th_Pmax:
            continue
        # Search for local maximum
        if wf[i] > wf[i-1] and wf[i] > wf[i+1]:
            return i, wf[i]
    return []

def find_second_peak(wf, th_Pmax_02):
    """
    Find the second peak above th_Pmax_02
    """

    for i in range(1, len(wf)):
        # Check on Pmax_02
        if wf[i] <= th_Pmax_02:
            continue
        # Search for local maximum
        if wf[i] > wf[i-1] and wf[i] > wf[i+1]:
            if wf[i] == max(wf):
                continue
            return i, wf[i]
    return []



def tfmra_retracker(P, th_rtck, th_P_01, th_P_02, gates_Pn):
    """
    Return the retracking gate and the index of the selected peak

    Input:
        - P             waveform power
        - th_rtck       Ppeak percentage to retrack
        - th_P          Pmax percentage to avoid retracking noise in front of
                        the leading edge
        - gates_Pn      number of gates used to compute thermal noise
    Output:
        - gate_rtck_01     retracking gate at peak 1
        - ipeak_01         index of the retracked peak 1
        - ipeakAmp_01      max amplitude of peak 1
        - th_tfmra01,      power treshold of peak 1
        - gate_rtck_02     retracking gate at peak 2
        - ipeak_02         index of the retracked peak 2
        - ipeakAmp_02      max amplitude of peak 2
        - th_tfmra02,      power treshold of peak 2
    """

    ### Pre-processing
    P       = list(P)
    gates   = range(len(P))
    AmpPmax_01    = max(P)                  # added by Henriette

    ### Retracking
    # Find first peak
    ipeak_01, ipeakAmp_01   = find_first_peak(P, max(P)*th_P_01)
    # Find second peak
    ipeak_02, ipeakAmp_02   = find_second_peak(P, max(P)*th_P_02)

    # Check if no peak is found
    if ipeak_01 == []:
        print("No peak found!")
        gate_rtck_01 = 'NaN'
        ipeak_01 = 'NaN'
        AmpPmax_01='NaN' 
        i_peakMax_01='NaN' 
        return float(gate_rtck_01), float(ipeak_01), float(AmpPmax_01), float(i_peakMax_01)  

    # If first peak is in the first or last gates_Pn gates, WF is discarded
    if ipeak_01 < gates_Pn or ipeak_01 > len(P) - gates_Pn:
        print("First peak coincides with thermal noise!")
        gate_rtck_01 = 'NaN'
        ipeak_01 = 'NaN'
        AmpPmax_01='NaN' 
        i_peakMax_01='NaN' 
        return float(gate_rtck_01), float(ipeak_01), float(AmpPmax_01), float(i_peakMax_01)    

    # Saving first peak power value
    Ppeak_01   = P[ipeak_01]
    Ppeak_02   = P[ipeak_02]

    # Compute thermal noise
    Pn      = np.mean(P[:gates_Pn])

    # Defining the power threshold
    th_tfmra01 = Pn + th_rtck * (Ppeak_01-Pn)

    # Find the first gate exceeding the power threshold:
    # searched after gates_Pn because the retracked point cannot be at the
    # gates used to compute thermal noise
    Ptemp   = next(x for x in P[gates_Pn:] if x > th_tfmra01)
    itemp   = P.index(Ptemp)

    # Calculate the exact retracked gate 01
    gate_rtck_01 = gates[itemp-1] + (th_tfmra01-P[itemp-1]) / (P[itemp]-P[itemp-1])


    #Find local minimum between the two peaks
    min_idx=np.argmin(P[ipeak_01:ipeak_02])


    # Defining the second power threshold
    th_tfmra02 = Pn + th_rtck * (Ppeak_02-Pn) 
    start_idx=ipeak_01+min_idx
    Ptemp   = next(x for x in P[start_idx:] if x > th_tfmra02)
    itemp   = P.index(Ptemp)

     # Calculate the exact retracked gate 02
    gate_rtck_02 = gates[itemp-1] + (th_tfmra02-P[itemp-1]) / (P[itemp]-P[itemp-1])



    return gate_rtck_01, ipeak_01, ipeakAmp_01, th_tfmra01, gate_rtck_02, ipeak_02, ipeakAmp_02, th_tfmra02
