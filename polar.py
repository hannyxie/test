import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
#读取文件，我这里用的是vcf文件
f = open('heading.raw.vcf', "r")
lines = f.read().split('\n')
#找到header
i = 0
for line in lines:
    if re.match("##", line, flags=0):
        i = i +1
    else:
        break
df = pd.DataFrame(lines[i+1:])
df = df[0].str.split('\t', expand=True)
df_index = lines[i][1:].split('\t')
df.columns=df_index
#对一些数值及类型进行转换
df = df.replace(['chr01','chr02','chr03','chr04','chr05','chr06','chr07','chr08','chr09','chr10','chr11','chr12'],[1,2,3,4,5,6,7,8,9,10,11,12])
df = df.replace('./.','0/0')
df = df.replace('.', 0)
df.dropna(inplace=True)
df['POS'] = df['POS'].astype('str').astype('int')
df['CHROM'] = df['CHROM'].astype('str').astype('int')
#仅仅保留我们需要的列就可以

#对变异类型进行区分
df = df[['CHROM', 'POS','REF', 'ALT']]
def get_type(ref,alt):
    if len(ref) == 1 and len(alt)==1:
        type = 'snp'
    else:
        type = 'indel'
    return type
df['type']=df.apply(lambda x: get_type(x.REF,x.ALT), axis = 1)

#接下来是进行滑窗，会生成一个新的df
chorms_max = df.groupby('CHROM')['POS'].max()
df_window = pd.DataFrame()
for i in range(1,13):
    max = chorms_max[i]
    window_start = 0
    while window_start < max:
        df_window = df_window.append({'chrom':int(i),'window_start':window_start},ignore_index=True)
        window_start = window_start+100*1000
df_window.eval('window_end = window_start+100*1000',inplace=True)
df_window.eval('window_central = window_start+50*1000',inplace=True)

#计算每个窗口内snp和indel的数量
def snp_count(chrom,window_start,window_end,df):
    snp_count= ''
    data_df = df[df['CHROM'] == chrom]
    data_df = data_df[(data_df['POS']<window_end) & (data_df['POS']> window_start)]
    data_df = data_df[data_df['type'] == 'snp']
    snp_count = data_df.shape[0]
    if snp_count != 0:
        snp_count = snp_count/1500
    return snp_count
def indel_count(chrom,window_start,window_end,df):
    indel_count= ''
    data_df = df[df['CHROM'] == chrom]
    data_df = data_df[(data_df['POS']<window_end) & (data_df['POS']> window_start)]
    data_df = data_df[data_df['type'] == 'indel']
    indel_count = data_df.shape[0]
    if indel_count != 0:
        indel_count = indel_count/1500
    return indel_count
df_window['snp_count']=df_window.apply(lambda x: snp_count(x.chrom,x.window_start,x.window_end,df),axis = 1)
df_window['indel_count']=df_window.apply(lambda x: indel_count(x.chrom,x.window_start,x.window_end,df),axis = 1)
print(chorms_max)

#开始绘图
fig = plt.figure(figsize=(10, 6))
plt.style.use('ggplot')
ax = plt.subplot(111, projection='polar')
#设置为顺时针
ax.set_theta_direction(-1)
#正上方为0度
ax.set_theta_zero_location('N')
#绘制柱状图，角度对应位置，半径对应高度
#先画最外面一圈,留下约15度的缺口，写图标，然后每个染色体中间空一格
chrom_length = [43270824,35937251,36413820,35502336,29958035,31248788,29697618,28443023,23012618,23206949,29020344,27531654]
i = 0
sum = 0
left = []
while i < 12:
    left.append(sum)
    sum += chrom_length[i]
    i += 1
n = 0
new_left = []
while n < 12:
    ax.barh(height=0.2,width=(2-(1/12)-(1/36))*np.pi*chrom_length[n]/sum,y=4,left = np.pi/12+(2-(1/12)-(1/36))*np.pi*left[n]/sum+np.pi/360*n)
    new_left.append(np.pi/12+(2-(1/12)-(1/36))*np.pi*left[n]/sum+np.pi/360*n)
    n += 1
#位置，高度，宽度,离圆心的距离
print(new_left)
j = 0
while j < 12:
    data_df = df_window[df_window['chrom'] == j+1]
    ax.bar(np.linspace(new_left[j],new_left[j]+(2-(1/12)-(1/36))*np.pi*data_df.shape[0]/3738,data_df.shape[0]),data_df['snp_count'],width=(2-(1/12)-(1/36))*np.pi/3738,bottom = 2.5,color='green')
    ax.bar(np.linspace(new_left[j],new_left[j]+(2-(1/12)-(1/36))*np.pi*data_df.shape[0]/3738,data_df.shape[0]),data_df['indel_count'],width=(2-(1/12)-(1/36))*np.pi/3738,bottom = 1.5,color='orange')
    j = j+1
ax.text(np.pi/24,4,'A')
ax.text(np.pi/24,3,'B')
ax.text(np.pi/24,2,'C')
new_left.append(2*np.pi)
k = 0

plt.tight_layout()
plt.axis('off')
plt.savefig('circus.png',dpi = 480)
plt.show()
