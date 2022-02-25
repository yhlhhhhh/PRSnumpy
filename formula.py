import numpy as np


'''
函数名称: or2prs
函数作用: 计算用户的多基因风险评分
参数 user_data: 用户基因型数据  dict{string 'rsid':string 'genotype'}
参数 minor: 次等位基因  numpy.array[string](1D)
参数 _or: 优势比  numpy.array[float](1D)
参数 maf: 次等位基因频率  numpy.array[float](1D)
参数 rsid: SNP的rsid  numpy.array[string](1D)

function name: or2prs
function: Calculate user's polygenic risk score
parameter user_data: user's genotype data  dict{string 'rsid':string 'genotype'}
parameter minor: minor alleles  numpy.array[string](1D)
parameter _or: odds ratio  numpy.array[float](1D)
parameter maf: minor alleles frequency  numpy.array[float](1D)
parameter rsid: SNP's rsid  numpy.array[string](1D)
'''
def or2prs(user_data, minor, _or, maf, rsid):
    alleles = []
    _2maf = maf * 2
    for i in range(minor.shape[0]):
        # 位点未检出时对应键值对不在字典中
        try:
            alleles.append(user_data[rsid[i]].count(minor[i]))
        except KeyError:
            alleles.append(_2maf[i])
    prs_value = np.prod(pow(_or, alleles))
    population_prs_value = np.sum(pow(_or, _2maf)) + 1
    return prs_value / population_prs_value


'''
函数名称: beta2prs
函数作用: 计算用户的多基因风险评分
参数 user_data: 用户基因型数据  dict{string 'rsid':string 'genotype'}
参数 minor: 次等位基因  numpy.array[string](1D)
参数 beta: 风险指数  numpy.array[float](1D)
参数 maf: 次等位基因频率  numpy.array[float](1D)
参数 rsid: SNP的rsid  numpy.array[string](1D)

function name: or2prs
function: Calculate user's polygenic risk score
parameter user_data: user's genotype data  dict{string 'rsid':string 'genotype'}
parameter minor: minor alleles  numpy.array[string](1D)
parameter beta: Beta coefficient  numpy.array[float](1D)
parameter maf: minor alleles frequency  numpy.array[float](1D)
parameter rsid: SNP's rsid  numpy.array[string](1D)
'''
def beta2prs(user_data, minor, beta, maf, rsid):
    alleles = []
    _2maf = maf * 2
    for i in range(minor.shape[0]):
        # 位点未检出时对应键值对不在字典中
        try:
            alleles.append(user_data[rsid[i]].count(minor[i]))
        except KeyError:
            alleles.append(_2maf[i])
    prs_value = np.dot(beta, alleles)
    population_prs_value = np.dot(beta, _2maf)
    return prs_value / population_prs_value