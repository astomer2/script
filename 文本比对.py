from pprint import pprint
import openpyxl
from openpyxl import load_workbook
import pandas as pd

Input = r"C:\Users\123\Desktop\jobwork\subject\PRL\dock\docking_result.xlsx"

# 打开工作簿
workbook = openpyxl.load_workbook(Input)

# 获取第二个工作表
ws1 = workbook.worksheets[1]
ws2 = workbook.worksheets[2]

# 获取第一列的值
text1 = [cell.value for cell in ws1['A']]
text2 = [cell.value for cell in ws2['A']]

def compare_text(text1, text2):
    set1 = set(text1)
    set2 = set(text2)

    unique_in_text1 = set1 - set2
    unique_in_text2 = set2 - set1

    return unique_in_text1, unique_in_text2
    

# 比较并输出结果
result1, result2 = compare_text(text1, text2)

print("第一列文本中独有的行：")
pprint(result1)

print("第二列文本中独有的行：")
pprint(result2)


# 手动输入两列文本
'''
text1 = """
GPCKAAF
PCKAAFP
CKAAFPR
KAAFPRW
AAFPRWF
AFPRWFY
FPRWFYD
PRWFYDT
RWFYDTK
WFYDTKT
FYDTKTG
YDTKTGC
DTKTGCS
TKTGCSG
KTGCSGF
TGCSGFI
GCSGFIY
CSGFIYG
SGFIYGG
GFIYGGC
FIYGGCD
IYGGCDG
YGGCDGN
GGCDGNR
GCDGNRN
CDGNRNN
DGNRNNE
GNRNNER
NRNNERT
GPCKAAFP
PCKAAFPR
CKAAFPRW
KAAFPRWF
AAFPRWFY
AFPRWFYD
FPRWFYDT
PRWFYDTK
RWFYDTKT
WFYDTKTG
FYDTKTGC
YDTKTGCS
DTKTGCSG
TKTGCSGF
KTGCSGFI
TGCSGFIY
GCSGFIYG
CSGFIYGG
SGFIYGGC
GFIYGGCD
FIYGGCDG
IYGGCDGN
YGGCDGNR
GGCDGNRN
GCDGNRNN
CDGNRNNE
DGNRNNER
GNRNNERT
GPCKAAFPR
PCKAAFPRW
CKAAFPRWF
KAAFPRWFY
AAFPRWFYD
AFPRWFYDT
FPRWFYDTK
PRWFYDTKT
RWFYDTKTG
WFYDTKTGC
FYDTKTGCS
YDTKTGCSG
DTKTGCSGF
TKTGCSGFI
KTGCSGFIY
TGCSGFIYG
GCSGFIYGG
CSGFIYGGC
SGFIYGGCD
GFIYGGCDG
FIYGGCDGN
IYGGCDGNR
YGGCDGNRN
GGCDGNRNN
GCDGNRNNE
CDGNRNNER
DGNRNNERT
GPCKAAFPRW
PCKAAFPRWF
CKAAFPRWFY
KAAFPRWFYD
AAFPRWFYDT
AFPRWFYDTK
FPRWFYDTKT
PRWFYDTKTG
RWFYDTKTGC
WFYDTKTGCS
FYDTKTGCSG
YDTKTGCSGF
DTKTGCSGFI
TKTGCSGFIY
KTGCSGFIYG
TGCSGFIYGG
GCSGFIYGGC
CSGFIYGGCD
SGFIYGGCDG
GFIYGGCDGN
FIYGGCDGNR
IYGGCDGNRN
YGGCDGNRNN
GGCDGNRNNE
GCDGNRNNER
CDGNRNNERT
GSICLEP
SICLEPK
ICLEPKV
CLEPKVV
LEPKVVG
EPKVVGP
PKVVGPC
KVVGPCT
VVGPCTA
VGPCTAY
GPCTAYF
PCTAYFR
CTAYFRR
TAYFRRF
AYFRRFY
YFRRFYF
FRRFYFN
RRFYFNS
RFYFNSE
FYFNSET
YFNSETG
FNSETGK
NSETGKC
SETGKCT
ETGKCTP
TGKCTPF
GKCTPFI
KCTPFIY
CTPFIYG
TPFIYGG
PFIYGGC
FIYGGCE
IYGGCEG
YGGCEGN
GGCEGNG
GCEGNGN
CEGNGNN
EGNGNNF
GNGNNFE
NGNNFET
GNNFETL
NNFETLR
NFETLRA
FETLRAC
ETLRACR
TLRACRA
LRACRAI
RACRAIC
ACRAICR
CRAICRA
DSTITIRGYVR
GPCTAYFR
PCTAYFRR
CTAYFRRF
TAYFRRFY
AYFRRFYF
YFRRFYFN
FRRFYFNS
RRFYFNSE
RFYFNSET
FYFNSETG
YFNSETGK
FNSETGKC
NSETGKCT
SETGKCTP
ETGKCTPF
TGKCTPFI
GKCTPFIY
KCTPFIYG
CTPFIYGG
TPFIYGGC
PFIYGGCE
FIYGGCEG
IYGGCEGN
YGGCEGNG
GGCEGNGN
GCEGNGNN
CEGNGNNF
EGNGNNFE
GNGNNFET
NGNNFETL
GNNFETLR
NNFETLRA
NFETLRAC
FETLRACR
ETLRACRA
TLRACRAI
LRACRAIC
RACRAICR
ACRAICRA
GSICLEPKV
SICLEPKVV
ICLEPKVVG
CLEPKVVGP
LEPKVVGPC
EPKVVGPCT
PKVVGPCTA
KVVGPCTAY
VVGPCTAYF
VGPCTAYFR
GPCTAYFRR
PCTAYFRRF
CTAYFRRFY
TAYFRRFYF
AYFRRFYFN
YFRRFYFNS
FRRFYFNSE
RRFYFNSET
RFYFNSETG
FYFNSETGK
YFNSETGKC
FNSETGKCT
NSETGKCTP
SETGKCTPF
ETGKCTPFI
TGKCTPFIY
GKCTPFIYG
KCTPFIYGG
CTPFIYGGC
TPFIYGGCE
PFIYGGCEG
FIYGGCEGN
IYGGCEGNG
YGGCEGNGN
GGCEGNGNN
GCEGNGNNF
CEGNGNNFE
EGNGNNFET
GNGNNFETL
NGNNFETLR
GNNFETLRA
NNFETLRAC
NFETLRACR
FETLRACRA
ETLRACRAI
TLRACRAIC
LRACRAICR
RACRAICRA
GSICLEPKVV
SICLEPKVVG
ICLEPKVVGP
CLEPKVVGPC
LEPKVVGPCT
EPKVVGPCTA
PKVVGPCTAY
KVVGPCTAYF
VVGPCTAYFR
VGPCTAYFRR
GPCTAYFRRF
PCTAYFRRFY
CTAYFRRFYF
TAYFRRFYFN
AYFRRFYFNS
YFRRFYFNSE
FRRFYFNSET
RRFYFNSETG
RFYFNSETGK
FYFNSETGKC
YFNSETGKCT
FNSETGKCTP
NSETGKCTPF
SETGKCTPFI
ETGKCTPFIY
TGKCTPFIYG
GKCTPFIYGG
KCTPFIYGGC
CTPFIYGGCE
TPFIYGGCEG
PFIYGGCEGN
FIYGGCEGNG
IYGGCEGNGN
YGGCEGNGNN
GGCEGNGNNF
GCEGNGNNFE
CEGNGNNFET
EGNGNNFETL
GNGNNFETLR
NGNNFETLRA
GNNFETLRAC
NNFETLRACR
NFETLRACRA
FETLRACRAI
ETLRACRAIC
TLRACRAICR
LRACRAICRA
FRRFYFD
RRFYFDS
RFYFDSE
FYFDSET
YFDSETG
FDSETGK
DSETGKC
NNFETLH
NFETLHA
FETLHAC
ETLHACR
TLHACRA
LHACRAI
HACRAIC
GSICLEPK
SICLEPKV
ICLEPKVV
CLEPKVVG
LEPKVVGP
EPKVVGPC
PKVVGPCT
KVVGPCTA
VVGPCTAY
VGPCTAYF
YFRRFYFD
FRRFYFDS
RRFYFDSE
RFYFDSET
FYFDSETG
YFDSETGK
FDSETGKC
DSETGKCT
GNNFETLH
NNFETLHA
NFETLHAC
FETLHACR
ETLHACRA
TLHACRAI
LHACRAIC
HACRAICR
AYFRRFYFD
YFRRFYFDS
FRRFYFDSE
RRFYFDSET
RFYFDSETG
FYFDSETGK
YFDSETGKC
FDSETGKCT
DSETGKCTP
NGNNFETLH
GNNFETLHA
NNFETLHAC
NFETLHACR
FETLHACRA
ETLHACRAI
TLHACRAIC
LHACRAICR
HACRAICRA
TAYFRRFYFD
AYFRRFYFDS
YFRRFYFDSE
FRRFYFDSET
RRFYFDSETG
RFYFDSETGK
FYFDSETGKC
YFDSETGKCT
FDSETGKCTP
DSETGKCTPF
GNGNNFETLH
NGNNFETLHA
GNNFETLHAC
NNFETLHACR
NFETLHACRA
FETLHACRAI
ETLHACRAIC
TLHACRAICR
LHACRAICRA
GGCEGNS
GCEGNSY
CEGNSYV
EGNSYVD
GNSYVDE
NSYVDEK
SYVDEKL
YVDEKLH
VDEKLHA
DEKLHAC
EKLHACR
KLHACRA
YGGCEGNS
GGCEGNSY
GCEGNSYV
CEGNSYVD
EGNSYVDE
GNSYVDEK
NSYVDEKL
SYVDEKLH
YVDEKLHA
VDEKLHAC
DEKLHACR
EKLHACRA
KLHACRAI
IYGGCEGNS
YGGCEGNSY
GGCEGNSYV
GCEGNSYVD
CEGNSYVDE
EGNSYVDEK
GNSYVDEKL
NSYVDEKLH
SYVDEKLHA
YVDEKLHAC
VDEKLHACR
DEKLHACRA
EKLHACRAI
KLHACRAIC
FIYGGCEGNS
IYGGCEGNSY
YGGCEGNSYV
GGCEGNSYVD
GCEGNSYVDE
CEGNSYVDEK
EGNSYVDEKL
GNSYVDEKLH
NSYVDEKLHA
SYVDEKLHAC
YVDEKLHACR
VDEKLHACRA
DEKLHACRAI
EKLHACRAIC
KLHACRAICR
"""

text2 = """
YYCRVRGGRC
YYCRVRGGR
YYCRVRGG
YYCRVRG
YSMEHFRWGK
YSMEHFRWG
YSMEHFRW
YSMEHFR
YPNGAEDES
YPNGAEDE
YPNGAED
YCRVRGGRCA
YCRVRGGRC
YCRVRGGR
YCRVRGG
WGKPVGKKRR
WGKPVGKKR
WGKPVGKK
WGKPVGK
VYPNGAEDES
VYPNGAEDE
VYPNGAED
VYPNGAE
VVRPRTPLSA
VVRPRTPLS
VVRPRTPL
VVRPRTP
VRPRTPLSAP
VRPRTPLSA
VRPRTPLS
VRPRTPL
VRGGRCAVLS
VRGGRCAVL
VRGGRCAV
VRGGRCA
VLSCLPKEEQ
VLSCLPKEE
VLSCLPKE
VLSCLPK
VKVYPNGAED
VKVYPNGAE
VKVYPNGA
VKVYPNG
VGKKRRPVKV
VGKKRRPVK
VGKKRRPV
VGKKRRP
VATRNSCKPP
VATRNSCKP
VATRNSCK
VATRNSC
TRNSCKPPAP
TRNSCKPPA
TRNSCKPP
TRNSCKP
TRGRKCCRRK
TRGRKCCRR
TRGRKCCR
TRGRKCC
TPLSAPCVAT
TPLSAPCVA
TPLSAPCV
TPLSAPC
TLQKYYCRVR
TLQKYYCRV
TLQKYYCR
TLQKYYC
SYSMEHFRWG
SYSMEHFRW
SYSMEHFR
SYSMEHF
STRGRKCCRR
STRGRKCCR
STRGRKCC
STRGRKC
SMEHFRWGKP
SMEHFRWGK
SMEHFRWG
SMEHFRW
SCRVLSLNC
SCRVLSLN
SCRVLSL
SCQCRFFRSA
SCQCRFFRS
SCQCRFFR
SCQCRFF
SCLPKEEQIG
SCLPKEEQI
SCLPKEEQ
SCLPKEE
SCKPPAPACC
SCKPPAPAC
SCKPPAPA
SCKPPAP
SAPCVATRNS
SAPCVATRN
SAPCVATR
SAPCVAT
SACSCRVLSL
SACSCRVLS
SACSCRVL
SACSCRV
RWGKPVGKKR
RWGKPVGKK
RWGKPVGK
RWGKPVG
RVRGGRCAVL
RVRGGRCAV
RVRGGRCA
RVRGGRC
RVLSLNC
RTPLSAPCVA
RTPLSAPCV
RTPLSAPC
RTPLSAP
RSACSCRVLS
RSACSCRVL
RSACSCRV
RSACSCR
RRPVKVYPNG
RRPVKVYPN
RRPVKVYP
RRPVKVY
RPVKVYPNGA
RPVKVYPNG
RPVKVYPN
RPVKVYP
RPRTPLSAPC
RPRTPLSAP
RPRTPLSA
RPRTPLS
RNSCKPPAPA
RNSCKPPAP
RNSCKPPA
RNSCKPP
RKCCRRKK
RKCCRRK
RGRKCCRRKK
RGRKCCRRK
RGRKCCRR
RGRKCCR
RGGRCAVLSC
RGGRCAVLS
RGGRCAVL
RGGRCAV
RFFRSACSCR
RFFRSACSC
RFFRSACS
RFFRSAC
RCAVLSCLPK
RCAVLSCLP
RCAVLSCL
RCAVLSC
QKYYCRVRGG
QKYYCRVRG
QKYYCRVR
QKYYCRV
QIGKCSTRGR
QIGKCSTRG
QIGKCSTR
QIGKCST
QCRFFRSACS
QCRFFRSAC
QCRFFRSA
QCRFFRS
PVKVYPNGAE
PVKVYPNGA
PVKVYPNG
PVKVYPN
PVGKKRRPVK
PVGKKRRPV
PVGKKRRP
PVGKKRR
PRTPLSAPCV
PRTPLSAPC
PRTPLSAP
PRTPLSA
PPAPACCDPC
PPAPACCDP
PPAPACCD
PPAPACC
PNGAEDES
PNGAEDE
PLSAPCVATR
PLSAPCVAT
PLSAPCVA
PLSAPCV
PKEEQIGKCS
PKEEQIGKC
PKEEQIGK
PKEEQIG
PCVATRNSCK
PCVATRNSC
PCVATRNS
PCVATRN
PCASCQCRFF
PCASCQCRF
PCASCQCR
PCASCQC
PAPACCDPCA
PAPACCDPC
PAPACCDP
PAPACCD
PACCDPCASC
PACCDPCAS
PACCDPCA
PACCDPC
NTLQKYYCRV
NTLQKYYCR
NTLQKYYC
NTLQKYY
NSCKPPAPAC
NSCKPPAPA
NSCKPPAP
NSCKPPA
NGAEDES
MEHFRWGKPV
MEHFRWGKP
MEHFRWGK
MEHFRWG
LSCLPKEEQI
LSCLPKEEQ
LSCLPKEE
LSCLPKE
LSAPCVATRN
LSAPCVATR
LSAPCVAT
LSAPCVA
LQKYYCRVRG
LQKYYCRVR
LQKYYCRV
LQKYYCR
LPKEEQIGKC
LPKEEQIGK
LPKEEQIG
LPKEEQI
KYYCRVRGGR
KYYCRVRGG
KYYCRVRG
KYYCRVR
KVYPNGAEDE
KVYPNGAED
KVYPNGAE
KVYPNGA
KVVRPRTPLS
KVVRPRTPL
KVVRPRTP
KVVRPRT
KRRPVKVYPN
KRRPVKVYP
KRRPVKVY
KRRPVKV
KPVGKKRRPV
KPVGKKRRP
KPVGKKRR
KPVGKKR
KPPAPACCDP
KPPAPACCD
KPPAPACC
KPPAPAC
KKVVRPRTPL
KKVVRPRTP
KKVVRPRT
KKVVRPR
KKRRPVKVYP
KKRRPVKVY
KKRRPVKV
KKRRPVK
KEEQIGKCST
KEEQIGKCS
KEEQIGKC
KEEQIGK
KCSTRGRKCC
KCSTRGRKC
KCSTRGRK
KCSTRGR
KCCRRKK
INTLQKYYCR
INTLQKYYC
INTLQKYY
INTLQKY
IINTLQKYYC
IINTLQKYY
IINTLQKY
IINTLQK
IGKCSTRGRK
IGKCSTRGR
IGKCSTRG
IGKCSTR
HFRWGKPVGK
HFRWGKPVG
HFRWGKPV
HFRWGKP
GRKCCRRKK
GRKCCRRK
GRKCCRR
GRCAVLSCLP
GRCAVLSCL
GRCAVLSC
GRCAVLS
GKPVGKKRRP
GKPVGKKRR
GKPVGKKR
GKPVGKK
GKKRRPVKVY
GKKRRPVKV
GKKRRPVK
GKKRRPV
GKCSTRGRKC
GKCSTRGRK
GKCSTRGR
GKCSTRG
GIINTLQKYY
GIINTLQKY
GIINTLQK
GIINTLQ
GGRCAVLSCL
GGRCAVLSC
GGRCAVLS
GGRCAVL
FRWGKPVGKK
FRWGKPVGK
FRWGKPVG
FRWGKPV
FRSACSCRVL
FRSACSCRV
FRSACSCR
FRSACSC
FFRSACSCRV
FFRSACSCR
FFRSACSC
FFRSACS
EQIGKCSTRG
EQIGKCSTR
EQIGKCST
EQIGKCS
EHFRWGKPVG
EHFRWGKPV
EHFRWGKP
EHFRWGK
EEQIGKCSTR
EEQIGKCST
EEQIGKCS
EEQIGKC
DPCASCQCRF
DPCASCQCR
DPCASCQC
DPCASCQ
CVATRNSCKP
CVATRNSCK
CVATRNSC
CVATRNS
CSTRGRKCCR
CSTRGRKCC
CSTRGRKC
CSTRGRK
CSCRVLSLNC
CSCRVLSLN
CSCRVLSL
CSCRVLS
CRVRGGRCAV
CRVRGGRCA
CRVRGGRC
CRVRGGR
CRVLSLNC
CRVLSLN
CRFFRSACSC
CRFFRSACS
CRFFRSAC
CRFFRSA
CQCRFFRSAC
CQCRFFRSA
CQCRFFRS
CQCRFFR
CLPKEEQIGK
CLPKEEQIG
CLPKEEQI
CLPKEEQ
CKPPAPACCD
CKPPAPACC
CKPPAPAC
CKPPAPA
CDPCASCQCR
CDPCASCQC
CDPCASCQ
CDPCASC
CCDPCASCQC
CCDPCASCQ
CCDPCASC
CCDPCAS
CAVLSCLPKE
CAVLSCLPK
CAVLSCLP
CAVLSCL
CASCQCRFFR
CASCQCRFF
CASCQCRF
CASCQCR
AVLSCLPKEE
AVLSCLPKE
AVLSCLPK
AVLSCLP
ATRNSCKPPA
ATRNSCKPP
ATRNSCKP
ATRNSCK
ASCQCRFFRS
ASCQCRFFR
ASCQCRFF
ASCQCRF
APCVATRNSC
APCVATRNS
APCVATRN
APCVATR
APACCDPCAS
APACCDPCA
APACCDPC
APACCDP
ACSCRVLSLN
ACSCRVLSL
ACSCRVLS
ACSCRVL
ACCDPCASCQ
ACCDPCASC
ACCDPCAS
ACCDPCA
"""
'''



