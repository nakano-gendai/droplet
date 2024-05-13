"""
Line plot etc. LaTeX font
10/01/2023
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import math
import matplotlib.ticker

### 入出力の設定
fsize = 0   # 0：小さいサイズ（論文用）　1：大きいサイズ（プレゼン用）
ext = 'pdf' # 保存ファイルの拡張子　pdf,svg,pngなどp
plotdir = 'plot2' # 保存用ディレクトリ
os.makedirs(plotdir, exist_ok = True) # ディレクトリ作成

# datadir = '/Users/nakanogendai/droplet/'
# datadir ='/home/nakano/anime/b4/data/data_fig6/'
datadir ="./"
plotdir, datadir, ext = plotdir + '/', datadir + '/', '.' + ext
dfile1   = datadir + 'dsd_test_3.d'
# dfile2   = datadir + 'taylortest2.d'
# dfile3   = datadir + 're0.2ca0.3.d'
# dfile4   = datadir + 're0.2ca0.4.d'
# dfile5   = datadir + 'ini_re0.2ca0.1.d'
# dfile6   = datadir + 'ini_re0.2ca0.2.d'
# dfile7   = datadir + 'ini_re0.2ca0.3.d'
# dfile8   = datadir + 'ini_re0.2ca0.4.d'



if os.path.exists(plotdir) == False or os.path.exists(datadir) == False: # plotdirやdatadirが存在しないときに中止する
    sys.exit("Error")

def sns_set(fs, tck_s, alw, ctxt):
    sns.set(
        context = ctxt,  # フォントサイズ・線幅（'paper' or 'talk'）
        palette = sns.color_palette("colorblind"),
        font = "serif",       # フォント指定
        # font = "sans-serif",  # サンセリフ体
        font_scale = fs,  # フォントスケール指定（これを変えるとcontextで決まるプリセットを更にいじれる）
        style = None, 
        # style = "whitegrid",   # 背景白，グリッドあり
        rc = {"text.usetex": True, 
        'text.latex.preamble' : r'\newcommand{\ssmr}[1]{_\mathrm{#1}}', # LaTeXのプリアンブル
        # 'text.latex.preamble' : r'\usepackage{txfonts}',   # nature系のフォント、ギリシャ文字あるとき
        'grid.linestyle': '--', 'grid.linewidth': 0, 
        "xtick.direction":"in", "xtick.major.width":0.8*alw, "xtick.major.size":tck_s, 
        "ytick.direction":"in", "ytick.major.width":0.8*alw, "ytick.major.size":tck_s, 
        "axes.linewidth":alw
        }
    )

### 読み込み
x1 = np.loadtxt(dfile1, usecols = 0, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
# x2 = np.loadtxt(dfile2, usecols = 0, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
# x3 = np.loadtxt(dfile3, usecols = 0, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
# x4 = np.loadtxt(dfile4, usecols = 0, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
# x5 = np.loadtxt(dfile5, usecols = 0, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
# x6 = np.loadtxt(dfile6, usecols = 0, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
# x7 = np.loadtxt(dfile7, usecols = 0, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
# x8 = np.loadtxt(dfile8, usecols = 0, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数

y1 = np.loadtxt(dfile1, usecols = 1, dtype = 'float64')
# y2 = np.loadtxt(dfile2, usecols = 1, dtype = 'float64')
# y3 = np.loadtxt(dfile3, usecols = 3, dtype = 'float64')
# y4 = np.loadtxt(dfile4, usecols = 3, dtype = 'float64')
# y5 = np.loadtxt(dfile5, usecols = 3, dtype = 'float64')
# y6 = np.loadtxt(dfile6, usecols = 3, dtype = 'float64')
# y7 = np.loadtxt(dfile7, usecols = 3, dtype = 'float64')
# y8 = np.loadtxt(dfile8, usecols = 3, dtype = 'float64')


### サイズ、ラベルなどの設定
# lx, ly = r'$step$', r'$\frac{\sum{V_{\mathrm{d}}}}{\sum{V_{\mathrm{d0}}}}$' # r''でTeX文字にできる
lx, ly = r'$D_{\mathrm{d}}/\eta_{\mathrm{k}}$', r'$P(D_{\mathrm{d}}/\eta_{\mathrm{k}})$' # r''でTeX文字にできる
if fsize == 0:
    fs1, lw1, ms1 = 1., 1., 2.8
    tck_s1, alw = 3, 0.625
    ctxt1 = 'paper'
    lpad = [5, 5] # 軸とラベルの間隔
    tpad = [5, 5] # 軸と数値の間隔
else:
    fs1, lw1, ms1 = 1.5, 3., 13.
    tck_s1, alw = 7, 1.25
    ctxt1, ext = 'talk', '_talk' + ext
    lpad = [10, 8]
    tpad = [12,20]


def plot_y1y2(): # y1, y2プロット用
    if fsize == 0:
        fig = plt.figure(figsize = (3.0, 3.0), dpi = 100, linewidth = 0)
    else:
        fig = plt.figure(figsize = (7.2, 7.0), dpi = 100, linewidth = 0)
    ax1 = fig.add_subplot(111)
    ax1.spines["top"].set_linewidth(alw)
    ax1.spines["left"].set_linewidth(alw)
    ax1.spines["bottom"].set_linewidth(alw)
    ax1.spines["right"].set_linewidth(alw)

    ax1.set_xlabel(lx, labelpad = lpad[0]) # 軸ラベル
    ax1.set_ylabel(ly, labelpad = lpad[1])
    xm, ym = [3, 250], [0.00001, 0.8]
    ax1.set_xlim(xm[0], xm[1]) # 軸の範囲
    ax1.set_ylim(ym[0], ym[1])
    # ax1.set_xticks(np.arange(0.0, xm[1] + 0.001, )) # xmaxまで0.2刻みの目盛り線
    ax1.tick_params(axis='x', pad = tpad[0])
    ax1.tick_params(axis='y', pad = tpad[1])
    plt.xscale('log')
    plt.yscale('log')
    ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    pos = [1, 10, 10**2] 
    ticks = ['1', '10', r'$10^2$']
    # pos = [10, 10**2] 
    # ticks = ['10', r'$10^2$']
    ax1.set_xticks(pos)
    ax1.set_xticklabels(ticks)

    ### 水平線と鉛直線
    ax1.axvline(25.85, lw = lw1*1.5, ls = 'dashdot', dashes = [1, 2], color = 'red', label = r'Kolmogorov--Hinze length $d_{\mathrm{H}}$') # dashesで破線の間隔などを設定できる

    za=np.arange(1.0, 60, 0.01)
    ua=1.7*10.0**(1.4)*za**(-3/2)
    ax1.plot(za, ua, lw = lw1*1.0, ls = 'dashed',  color = 'black', 
            alpha = 1,     # 透明度
            clip_on = True, # プロット枠外にもプロットする
            zorder = 9,      # zorderが大きいほど前面に表示される
            label = r'$\propto D_{\mathrm{d}}^{-3/2}$') 
    
    za=np.arange(35, 250, 0.01)
    ua=8000.0**(1.15)*za**(-10/3)
    ax1.plot(za, ua, lw = lw1*1.0, ls = 'dotted',  color = 'black', 
            alpha = 1,     # 透明度
            clip_on = True, # プロット枠外にもプロットする
            zorder = 9,      # zorderが大きいほど前面に表示される
            label = r'$\propto D_{\mathrm{d}}^{-10/3}$') 
    
    # za=np.arange(1.0, 60, 0.01)
    # ua=1.7*10.0**(1.15)*za**(-3/2)
    # ax1.plot(za, ua, lw = lw1*1.0, ls = 'dashed',  color = 'black', 
    #         alpha = 1,     # 透明度
    #         clip_on = True, # プロット枠外にもプロットする
    #         zorder = 9,      # zorderが大きいほど前面に表示される
    #         label = r'$\propto D_{\mathrm{d}}^{-3/2}$') 
    
    # za=np.arange(30, 250, 0.01)
    # ua=8000.0**(1.08)*za**(-10/3)
    # ax1.plot(za, ua, lw = lw1*1.0, ls = 'dotted',  color = 'black', 
    #         alpha = 1,     # 透明度
    #         clip_on = True, # プロット枠外にもプロットする
    #         zorder = 9,      # zorderが大きいほど前面に表示される
    #         label = r'$\propto D_{\mathrm{d}}^{-10/3}$') 
    
    # ax1.plot(x1, y1, lw = 1, ls = 'none', marker = 'o', ms = ms1*2.5, mew = lw1*1.5, mfc = 'none', color = 'blue', 
    #         alpha = 1,     # 透明度
    #         clip_on = False, # プロット枠外にもプロットする
    #         zorder = 1,      # zorderが大きいほど前面に表示される
    #         label = 'Probability density function') 
    ax1.plot(x1, y1, lw = lw1*2.5, ls = 'solid',  color = 'gray', alpha = 0.7, clip_on = True, zorder = 13, label = 'Probability density function')

    # 凡例の設定
    h1, l1 = ax1.get_legend_handles_labels()
    ax1.legend(h1, l1, 
    bbox_to_anchor = (1.0, 1.0), loc = "upper left", # bbox_to_anchorは凡例のlocの座標
    framealpha = 1.0, fancybox=False, fontsize=8.0,
    edgecolor = "black").get_frame().set_linewidth(alw*0.8)
    ## 保存
    fig.savefig(plotdir + "DSD_test_2_3" + ext, bbox_inches = "tight") # bbox_inches="tight"で余白をなくす

##=================== main ===================##
if __name__ == '__main__':
    print("start main")
    print("datadir: ", datadir, ", plotdir: ", plotdir, ", ctxt: ", ctxt1)
    sns_set(fs1, tck_s1, alw, ctxt1)
    plot_y1y2()
    print("end main")