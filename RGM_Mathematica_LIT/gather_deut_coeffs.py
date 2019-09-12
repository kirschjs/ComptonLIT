import os, re
from math import ceil


def parse_ev_coeffs(mult=0, outf='deut-coeffs-lit'):
    out = [line2 for line2 in open('OUTPUT')]
    #for n in range(1,len(out)):
    #    if(out[n].strip()=="EIGENWERTE DES HAMILTONOPERATORS"):
    #        print(float(out[n+3].split()[0]))
    coef = ''
    coeffp = []
    coeff_mult = []
    bvc = 0
    for line in range(0, len(out) - 1):
        if re.search('ENTWICKLUNG DES  1 TEN EIGENVEKTORS', out[line]):
            for bvl in range(line + 2, len(out)):
                if ((out[bvl][:3] == ' KO') | (out[bvl][:3] == '\n')):
                    bvc = out[bvl - 1].strip().split('/')[-1].split(')')[0]
                    break
                coeffp += [
                    float(coo.split('/')[0])
                    for coo in out[bvl].strip().split(')')[:-1]
                ]
                coef += out[bvl]
            break
    s = ''
    for n in range(len(coeffp)):
        if mult:
            for m in range(len(coeffp) - n):
                if m == 0:
                    s += '%18.10g' % (coeffp[n] * coeffp[n + m]) + '\n'
                # for identical fragments, c1*c2|BV1>|BV2> appears twice and can be summed up => faktor 2
                # see coef_mul_id.exe
                else:
                    s += '%18.10g' % (coeffp[n] * coeffp[n + m] * 2) + '\n'
        else:
            s += '%E' % (coeffp[n]) + '\n'
            #s += '%18.10g' % (coeffp[n]) + '\n'
    ss = s.replace('e', 'E')
    if bvc == 0:
        print("No coefficients found in OUTPUT")
    with open(outf, 'w') as outfile:
        outfile.write(ss)
    return


def parse_widths(outf='gauss-widths-lit'):
    out = [line2 for line2 in open('INQUA_N')]

    s = ''
    for wl in range(5, 5 + ceil(int(out[3].split()[1]) / 6)):
        print(out[wl])

        wline = out[wl].strip().split()
        for n in range(len(wline)):
            s += '%18.10g' % float(wline[n]) + '\n'

    with open(outf, 'w') as outfile:
        outfile.write(s)
    return


parse_ev_coeffs()
#parse_widths()
