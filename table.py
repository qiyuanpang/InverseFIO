
def oneD(n):
    filename1 = "results" + str(n) + "_1D_6.log"
    filename2 = "results" + str(n) + "_1D_3.log"
    data1 = []
    data2 = []
    f1 = open(filename1, mode='r')
    f2 = open(filename2, mode='r')
    for line in f1:
        line = line.split()
        data1.append(line)
    f1.close()
    for line in f2:
        line = line.split()
        data2.append(line)

    if n != 4:
        lines = []
        lines.append(['256', '1e-7', '1e-6', data1[53][2][0:10], data1[54][11], data1[54][15], data1[55][8], data1[55][12],
                                    '1e-3', data2[53][2][0:10], data2[54][11], data2[54][15], data2[55][8], data2[55][12]])
        lines.append(['625', '1e-7', '1e-6', data1[64][2][0:10], data1[65][11], data1[65][15], data1[66][8], data1[66][12],
                                    '1e-3', data2[64][2][0:10], data2[65][11], data2[65][15], data2[66][8], data2[66][12]])
        lines.append(['1296', '1e-7','1e-6', data1[75][2][0:10], data1[76][11], data1[76][15], data1[77][8], data1[77][12],
                                    '1e-3', data2[75][2][0:10], data2[76][11], data2[76][15], data2[77][8], data2[77][12]])
        lines.append(['2401', '1e-7','1e-6', data1[86][2][0:10], data1[87][11], data1[87][15], data1[88][8], data1[88][12],
                                    '1e-3', data2[86][2][0:10], data2[87][11], data2[87][15], data2[88][8], data2[88][12]])
        lines.append(['4096', '1e-7','1e-6', data1[97][2][0:10], data1[98][11], data1[98][15], data1[99][8], data1[99][12],
                                    '1e-3', data2[97][2][0:10], data2[98][11], data2[98][15], data2[99][8], data2[99][12]])
        lines.append(['6561', '1e-7','1e-6', data1[108][2][0:10], data1[109][11], data1[109][15], data1[110][8], data1[110][12],
                                    '1e-3', data2[108][2][0:10], data2[109][11], data2[109][15], data2[110][8], data2[110][12]])
    else:
        lines = []
        lines.append(['256', '1e-7', '1e-6', data1[53][2][0:10], data1[54][11], data1[54][15], data1[55][8], data1[55][12],
                                    '1e-3', data2[53][2][0:10], data2[54][11], data2[54][15], data2[55][8], data2[55][12]])
        lines.append(['625', '1e-7', '1e-6', data1[64][2][0:10], data1[65][11], data1[65][15], data1[66][8], data1[66][12],
                                    '1e-3', data2[64][2][0:10], data2[65][11], data2[65][15], data2[66][8], data2[66][12]])
        lines.append(['1296', '1e-7','1e-6', data1[75][2][0:10], data1[76][11], data1[76][15], data1[77][8], data1[77][12],
                                    '1e-3', data2[75][2][0:10], data2[76][11], data2[76][15], data2[77][8], data2[77][12]])
        lines.append(['2401', '1e-7','1e-6', data1[86][2][0:10], data1[87][11], data1[87][15], data1[88][8], data1[88][12],
                                    '1e-3', data2[86][2][0:10], data2[87][11], data2[87][15], data2[88][8], data2[88][12]])
        lines.append(['4096', '1e-7','1e-6', data1[97][2][0:10], data1[98][11], data1[98][15], data1[99][8], data1[99][12],
                                    '1e-3', data2[97][2][0:10], data2[98][11], data2[98][15], data2[99][8], data2[99][12]])
        lines.append(['6561', '1e-7','1e-6', data1[149][2][0:10], data1[650][11], data1[650][15], data1[651][8], data1[651][12],
                                    '1e-3', data2[149][2][0:10], data2[650][11], data2[650][15], data2[651][8], data2[651][12]])
    for line in lines:
        print('\hline')
        print('\multirow{2}*{'+line[0]+'} & \multirow{2}*{'+line[1]+'} & '+line[2]+' & '+line[3]+' & '+line[4]+' & '+line[5]+' & '+line[6]+' & '+line[7]+' \\\\')
        print('~ & ~ & '+line[8]+' & '+line[9]+' & '+line[10]+' & '+line[11]+' & '+line[12]+' & '+line[13]+' \\\\')
                            
def oneDf(n):
    filename1 = "results" + str(n) + "_1D_f6.log"
    filename2 = "results" + str(n) + "_1D_f3.log"
    data1 = []
    data2 = []
    f1 = open(filename1, mode='r')
    f2 = open(filename2, mode='r')
    for line in f1:
        line = line.split()
        data1.append(line)
    f1.close()
    for line in f2:
        line = line.split()
        data2.append(line)
    f2.close()
    if n != 2:
        lines = []
        lines.append(['256', '1e-7', '$15 log_{2}N$', data1[53][2][0:10], data1[54][11], data1[54][15], data1[55][8], data1[55][12],
                                    '$8 log_{2}N$', data2[53][2][0:10], data2[54][11], data2[54][15], data2[55][8], data2[55][12]])
        lines.append(['625', '1e-7', '$15 log_{2}N$', data1[69][2][0:10], data1[70][11], data1[70][15], data1[71][8], data1[71][12],
                                    '$8 log_{2}N$', data2[69][2][0:10], data2[70][11], data2[70][15], data2[71][8], data2[71][12]])
        lines.append(['1296', '1e-7','$15 log_{2}N$', data1[85][2][0:10], data1[86][11], data1[86][15], data1[87][8], data1[87][12],
                                    '$8 log_{2}N$', data2[85][2][0:10], data2[86][11], data2[86][15], data2[87][8], data2[87][12]])
        lines.append(['2401', '1e-7','$15 log_{2}N$', data1[101][2][0:10], data1[102][11], data1[102][15], data1[103][8], data1[103][12],
                                    '$8 log_{2}N$', data2[101][2][0:10], data2[102][11], data2[102][15], data2[103][8], data2[103][12]])
        lines.append(['4096', '1e-7','$15 log_{2}N$', data1[117][2][0:10], data1[118][11], data1[118][15], data1[119][8], data1[119][12],
                                    '$8 log_{2}N$', data2[117][2][0:10], data2[118][11], data2[118][15], data2[119][8], data2[119][12]])
        lines.append(['6561', '1e-7','$15 log_{2}N$', data1[133][2][0:10], data1[134][11], data1[134][15], data1[135][8], data1[135][12],
                                    '$8 log_{2}N$', data2[133][2][0:10], data2[134][11], data2[134][15], data2[135][8], data2[135][12]])
    else:
        lines = []
        lines.append(['256', '1e-7', '$15 log_{2}N$', data1[53][2][0:10], data1[54][11], data1[54][15], data1[55][8], data1[55][12],
                                    '$8 log_{2}N$', data2[53][2][0:10], data2[58][11], data2[58][15], data2[59][8], data2[59][12]])
        lines.append(['625', '1e-7', '$15 log_{2}N$', data1[69][2][0:10], data1[70][11], data1[70][15], data1[71][8], data1[71][12],
                                    '$8 log_{2}N$', data2[73][2][0:10], data2[74][11], data2[74][15], data2[75][8], data2[75][12]])
        lines.append(['1296', '1e-7','$15 log_{2}N$', data1[85][2][0:10], data1[86][11], data1[86][15], data1[87][8], data1[87][12],
                                    '$8 log_{2}N$', data2[89][2][0:10], data2[90][11], data2[90][15], data2[91][8], data2[91][12]])
        lines.append(['2401', '1e-7','$15 log_{2}N$', data1[101][2][0:10], data1[102][11], data1[102][15], data1[103][8], data1[103][12],
                                    '$8 log_{2}N$', data2[105][2][0:10], data2[106][11], data2[106][15], data2[107][8], data2[107][12]])
        lines.append(['4096', '1e-7','$15 log_{2}N$', data1[117][2][0:10], data1[118][11], data1[118][15], data1[119][8], data1[119][12],
                                    '$8 log_{2}N$', data2[121][2][0:10], data2[122][11], data2[122][15], data2[123][8], data2[123][12]])
        lines.append(['6561', '1e-7','$15 log_{2}N$', data1[133][2][0:10], data1[134][11], data1[134][15], data1[135][8], data1[135][12],
                                    '$8 log_{2}N$', data2[137][2][0:10], data2[138][11], data2[138][15], data2[139][8], data2[139][12]])


    for line in lines:
        print('\hline')
        print('\multirow{2}*{'+line[0]+'} & \multirow{2}*{'+line[1]+'} & '+line[2]+' & '+line[3]+' & '+line[4]+' & '+line[5]+' & '+line[6]+' & '+line[7]+' \\\\')
        print('~ & ~ & '+line[8]+' & '+line[9]+' & '+line[10]+' & '+line[11]+' & '+line[12]+' & '+line[13]+' \\\\')
  
def twoD(n):
    filename1 = "results" + str(n) + "_2D_3.log"
    data1 = []
    f1 = open(filename1, mode='r')
    for line in f1:
        line = line.split()
        data1.append(line)
    f1.close()
    if n != 4:
        lines = []
        lines.append(['$16^2$', '1e-7', '1e-3', data1[48][2][0:10], data1[49][11], data1[49][15], data1[50][8], data1[50][12]])
        lines.append(['$25^2$', '1e-7', '1e-3', data1[59][2][0:10], data1[60][11], data1[60][15], data1[61][8], data1[61][12]])
        lines.append(['$36^2$', '1e-7', '1e-3', data1[70][2][0:10], data1[71][11], data1[71][15], data1[72][8], data1[72][12]])
        lines.append(['$49^2$', '1e-7', '1e-3', data1[81][2][0:10], data1[82][11], data1[82][15], data1[83][8], data1[83][12]])
        lines.append(['$64^2$', '1e-7', '1e-3', data1[92][2][0:10], data1[93][11], data1[93][15], data1[94][8], data1[94][12]])
        lines.append(['$81^2$', '1e-7', '1e-3', data1[103][2][0:10], data1[104][11], data1[104][15], data1[105][8], data1[105][12]])
    else:
        lines = []
        lines.append(['$16^2$', '1e-7', '1e-3', data1[48][2][0:10], data1[49][11], data1[49][15], data1[50][8], data1[50][12]])
        lines.append(['$25^2$', '1e-7', '1e-3', data1[59][2][0:10], data1[60][11], data1[60][15], data1[61][8], data1[61][12]])
        lines.append(['$36^2$', '1e-7', '1e-3', data1[70][2][0:10], data1[75][11], data1[75][15], data1[76][8], data1[76][12]])
        lines.append(['$49^2$', '1e-7', '1e-3', data1[85][2][0:10], data1[86][11], data1[86][15], data1[87][8], data1[87][12]])
        lines.append(['$64^2$', '1e-7', '1e-3', data1[96][2][0:10], data1[97][11], data1[97][15], data1[98][8], data1[98][12]])
        lines.append(['$81^2$', '1e-7', '1e-3', data1[107][2][0:10], data1[108][11], data1[108][15], data1[109][8], data1[109][12]])

    for line in lines:
        print('\hline')
        print(line[0]+' & '+line[1]+' & '+line[2]+' & '+line[3]+' & '+line[4]+' & '+line[5]+' & '+line[6]+' & '+line[7]+' \\\\')

def twoDf(n):
    filename1 = "results" + str(n) + "_2D.log"
    data1 = []
    f1 = open(filename1, mode='r')
    for line in f1:
        line = line.split()
        data1.append(line)
    f1.close()
    if n != 4:
        lines = []
        lines.append(['$16^2$', '1e-7', '$8 log_{2}N$', data1[48][2][0:10], data1[49][11], data1[49][15], data1[50][8], data1[50][12]])
        lines.append(['$25^2$', '1e-7', '$8 log_{2}N$', data1[64][2][0:10], data1[65][11], data1[65][15], data1[66][8], data1[66][12]])
        lines.append(['$36^2$', '1e-7', '$8 log_{2}N$', data1[80][2][0:10], data1[81][11], data1[81][15], data1[82][8], data1[82][12]])
        lines.append(['$49^2$', '1e-7', '$8 log_{2}N$', data1[96][2][0:10], data1[97][11], data1[97][15], data1[98][8], data1[98][12]])
        lines.append(['$64^2$', '1e-7', '$8 log_{2}N$', data1[112][2][0:10], data1[113][11], data1[113][15], data1[114][8], data1[114][12]])
        lines.append(['$81^2$', '1e-7', '$8 log_{2}N$', data1[128][2][0:10], data1[129][11], data1[129][15], data1[130][8], data1[130][12]])
    else:
        lines = []
        lines.append(['$16^2$', '1e-7', '$8 log_{2}N$', data1[48][2][0:10], data1[49][11], data1[49][15], data1[50][8], data1[50][12]])
        lines.append(['$25^2$', '1e-7', '$8 log_{2}N$', data1[64][2][0:10], data1[65][11], data1[65][15], data1[66][8], data1[66][12]])
        lines.append(['$36^2$', '1e-7', '$8 log_{2}N$', data1[80][2][0:10], data1[81][11], data1[81][15], data1[82][8], data1[82][12]])
        lines.append(['$49^2$', '1e-7', '$8 log_{2}N$', data1[96][2][0:10], data1[97][11], data1[97][15], data1[98][8], data1[98][12]])
        lines.append(['$64^2$', '1e-7', '$8 log_{2}N$', data1[112][2][0:10], data1[113][11], data1[113][15], data1[114][8], data1[114][12]])
        lines.append(['$81^2$', '1e-7', '$8 log_{2}N$', data1[128][2][0:10], data1[129][11], data1[129][15], data1[130][8], data1[130][12]])

    for line in lines:
        print('\hline')
        print(line[0]+' & '+line[1]+' & '+line[2]+' & '+line[3]+' & '+line[4]+' & '+line[5]+' & '+line[6]+' & '+line[7]+' \\\\')
   
 

if __name__ == '__main__':
    twoDf(4)
