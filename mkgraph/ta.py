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

### 入出力の設定
fsize = 0   # 0：小さいサイズ（論文用）　1：大きいサイズ（プレゼン用）
ext = 'pdf' # 保存ファイルの拡張子　pdf,svg,pngなどp
plotdir = './' # 保存用ディレクトリ
os.makedirs(plotdir, exist_ok = True) # ディレクトリ作成

# datadir = '/Users/nakanogendai/droplet/'
# datadir ='/home/nakano/anime/b4/data/data_fig6/'
datadir ="./"
plotdir, datadir, ext = plotdir + '/', datadir + '/', '.' + ext
dfile1   = datadir + '30.d'
dfile2   = datadir + '40.d'
dfile3   = datadir + '50.d'
dfile4   = datadir + '60.d'
dfile5   = datadir + '70.d'
# dfile6   = datadir + 'case1_e.d'
# dfile7   = datadir + 'case2_e.d'
# dfile8   = datadir + 'case20_e.d'
# dfile9   = datadir + 'case21_e.d'
# dfile10   = datadir + 'case26_e.d'

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
x2 = np.loadtxt(dfile2, usecols = 0, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
x3 = np.loadtxt(dfile3, usecols = 0, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
x4 = np.loadtxt(dfile4, usecols = 0, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
x5 = np.loadtxt(dfile5, usecols = 0, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
# x6 = np.loadtxt(dfile1, usecols = 0, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
# x7 = np.loadtxt(dfile1, usecols = 0, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
# x8 = np.loadtxt(dfile1, usecols = 0, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
# x9 = np.loadtxt(dfile9, usecols = 0, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
# x10 = np.loadtxt(dfile10, usecols = 0, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数

y1 = np.loadtxt(dfile1, usecols = 1, dtype = 'float64')
y2 = np.loadtxt(dfile2, usecols = 1, dtype = 'float64')
y3 = np.loadtxt(dfile3, usecols = 1, dtype = 'float64')
y4 = np.loadtxt(dfile4, usecols = 1, dtype = 'float64')
y5 = np.loadtxt(dfile5, usecols = 1, dtype = 'float64')
# y6 = np.loadtxt(dfile1, usecols = 6, dtype = 'float64')
# y7 = np.loadtxt(dfile1, usecols = 7, dtype = 'float64')
# y8 = np.loadtxt(dfile1, usecols = 8, dtype = 'float64')
# y9 = np.loadtxt(dfile9, usecols = 2, dtype = 'float64')
# y10 = np.loadtxt(dfile10, usecols = 2, dtype = 'float64')



### サイズ、ラベルなどの設定
lx, ly = r'$\mathrm{We}$', r'$N_\mathrm{d}$' # r''でTeX文字にできる
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
    
    # xm, ym = [0, 2.8], [-0.000000000001, 0.00000000001]
    xm, ym = [0.4, 1.6], [0, 3]
    ax1.set_xlim(xm[0], xm[1]) # 軸の範囲
    ax1.set_ylim(ym[0], ym[1])
    # ax1.set_yticks(np.arange(0.0, ym[1] + 0.0000005)) # xmaxまで0.2刻みの目盛り線
    ax1.tick_params(axis='x', pad = tpad[0])
    ax1.tick_params(axis='y', pad = tpad[1])
    # plt.yscale('log')
    # plt.xscale('log')

    pos = [0.5, 1, 1.5] 
    ticks = [r'$0.5$', r'$1$', r'$1.5$']
    ax1.set_xticks(pos)
    ax1.set_xticklabels(ticks)
    
    pos = [0, 1, 2, 3] 
    ticks = [r'$0$', r'$1$', r'$2$', r'$3$']
    ax1.set_yticks(pos)
    ax1.set_yticklabels(ticks)
    
    # pos = [-0.25, 0, 20] 
    # ticks = [r'$0$', r'$10$', r'$20$']
    # ax1.set_yticks(pos)
    # ax1.set_yticklabels(ticks)
    
    # pos = [0.01, 0.1, 1] 
    # ticks = [r'$10^{-2}$', r'$10^{-1}$', r'$1$']
    # ax1.set_xticks(pos)
    # ax1.set_xticklabels(ticks)
    
    # pos = [1, 10, 100] 
    # ticks = [r'$1$', r'$10$', r'$10^{2}$']
    # ax1.set_xticks(pos)
    # ax1.set_xticklabels(ticks)
    
    # pos = [0, 5, 10, 15, 20] 
    # ticks = [r'$0$', r'$5$', r'$10$', r'$15$', r'$20$']
    # ax1.set_xticks(pos)
    # ax1.set_xticklabels(ticks)

    ### 水平線と鉛直線
    # ax1.axhline(0, lw = lw1*0.8, ls = 'dotted', label = r'$y=0$')
    # ax1.axvline(7.0, lw = lw1*0.8, ls = 'dashed', dashes = [1, 2], color = 'red', label = r'filter line') # dashesで破線の間隔などを設定できる
    # ytext = -22.
    # ax1.axhline(ytext, lw = lw1*0.8, c = 'C4', ls = 'dashdot')
    # ax1.text(0.7, (ytext - ym[0])/(ym[1] - ym[0]) + 0.05, r'$y=%4.1f$'%(ytext), transform = ax1.transAxes) # 図中にtextを入れる
    

    # za=np.arange(1, 260, 0.0001)
    # ua=2.5**(-10.5)*za**(-5/3)
    # ax1.plot(za, ua, lw = lw1*0.8, ls = 'dotted',  color = 'black', 
    #         alpha = 0.5,     # 透明度
    #         clip_on = True, # プロット枠外にもプロットする
    #         zorder = 9,      # zorderが大きいほど前面に表示される
    #         label = r'$\propto k^{-5/3}$') 

    
    ###color map
    # a = x1.reshape(66,150)
    # ax1.set_aspect("equal")
    # extent1=[-0.5+xm[0], 0.5+xm[1], -0.5+ym[0], 0.5+ym[1]]
    # im = ax1.imshow(a, origin="lower", cmap=plt.cm.plasma, vmin=2, vmax=5, extent=extent1)
    # axpos = ax1.get_position()
    # pp_ax = fig.add_axes([axpos.x1 + 0.02, axpos.y0, 0.05, axpos.height]) # 左下のx,y,幅,高さ?
    # pp=fig.colorbar(im, ax=ax1, orientation="vertical", cax=pp_ax)
    
    ### 読み込みデータのプロット
    ax1.plot(x1, y1, lw = 1, ls = 'none', marker = 'o', ms = ms1*3.0, mew = lw1*1, mfc = 'none', color = 'red', 
            alpha = 0.8,     # 透明度
            clip_on = False, # プロット枠外にもプロットする
            zorder = 110,
            label = r'$D/\eta_\mathrm{k}=30$') 
    ax1.plot(x2, y2, lw = 1, ls = 'none', marker = '^', ms = ms1*3.0, mew = lw1*1, mfc = 'none', color = 'black', 
            alpha = 0.8,     # 透明度
            clip_on = False, # プロット枠外にもプロットする
            zorder = 109,      # zorderが大きいほど前面に表示される
            label = r'$D/\eta_\mathrm{k}=40$') 
    ax1.plot(x3, y3, lw = 1, ls = 'none', marker = 'D', ms = ms1*3.0, mew = lw1*1, mfc = 'none', color = 'blue', 
            alpha = 0.8,     # 透明度
            clip_on = False, # プロット枠外にもプロットする
            zorder = 108,      # zorderが大きいほど前面に表示される
            label = r'$D/\eta_\mathrm{k}=50$') 
    ax1.plot(x4, y4, lw = 1, ls = 'none', marker = '*', ms = ms1*3.0, mew = lw1*1, mfc = 'none', color = 'gray', 
            alpha = 0.8,     # 透明度
            clip_on = False, # プロット枠外にもプロットする
            zorder = 107,      # zorderが大きいほど前面に表示される
            label = r'$D/\eta_\mathrm{k}=60$') 
    ax1.plot(x5, y5, lw = 1, ls = 'none', marker = '+', ms = ms1*3.0, mew = lw1*1, mfc = 'none', color = 'deepskyblue', 
            alpha = 0.8,     # 透明度
            clip_on = False, # プロット枠外にもプロットする
            zorder = 106,      # zorderが大きいほど前面に表示される
            label = r'$D/\eta_\mathrm{k}=70$') 
    # ax1.plot(x6, y6, lw = 1, ls = 'none', marker = 's', ms = ms1*1.0, mew = lw1*0.5, mfc = 'none', color = 'orange', 
    #         alpha = 0.8,     # 透明度
    #         clip_on = False, # プロット枠外にもプロットする
    #         zorder = 105,      # zorderが大きいほど前面に表示される
    #         label = r'$k=6$') 
    # ax1.plot(x7, y7, lw = 1, ls = 'none', marker = 'p', ms = ms1*1.0, mew = lw1*0.5, mfc = 'none', color = 'purple', 
    #         alpha = 0.8,     # 透明度
    #         clip_on = False, # プロット枠外にもプロットする
    #         zorder = 104,      # zorderが大きいほど前面に表示される
    #         label = r'$k=7$') 
    # ax1.plot(x8, y8, lw = 1, ls = 'none', marker = 'x', ms = ms1*1.0, mew = lw1*0.5, mfc = 'none', color = 'lime', 
    #         alpha = 0.8,     # 透明度
    #         clip_on = False, # プロット枠外にもプロットする
    #         zorder = 103,      # zorderが大きいほど前面に表示される
    #         label = r'$k=8$') 
    # ax1.plot(x8, y8, lw = 1, ls = 'none', marker = 'x', ms = ms1*1.0, mew = lw1*0.5, mfc = 'none', color = 'lime', 
    #         alpha = 0.8,     # 透明度
    #         clip_on = False, # プロット枠外にもプロットする
    #         zorder = 103,      # zorderが大きいほど前面に表示される
    #         label = r'$k=8$') 

    
    # ax1.plot(x1, y1, lw = lw1*0.5, ls = 'solid',  color = 'red', alpha = 0.8, clip_on = True, zorder = 20)
    # ax1.plot(x2, y2, lw = lw1*0.5, ls = 'solid',  color = 'black', alpha = 0.8, clip_on = True, zorder = 19)
    # ax1.plot(x3, y3,lw = lw1*0.5, ls = 'solid',  color = 'blue', alpha = 0.8, clip_on = True, zorder = 18)
    # ax1.plot(x4, y4, lw = lw1*0.5, ls = 'solid',  color = 'gray', alpha = 0.8, clip_on = True, zorder = 17)
    # ax1.plot(x5, y5, lw = lw1*0.5, ls = 'solid',  color = 'deepskyblue', alpha = 0.8, clip_on = True, zorder = 16)
    # ax1.plot(x6, y6, lw = lw1*0.5, ls = 'solid',  color = 'orange', alpha = 0.8, clip_on = True, zorder = 15)
    # ax1.plot(x7, y7, lw = lw1*0.5, ls = 'solid',  color = 'purple', alpha = 0.8, clip_on = True, zorder = 14)
    # ax1.plot(x8, y8, lw = lw1*0.5, ls = 'solid',  color = 'lime', alpha = 0.8, clip_on = True, zorder = 13)
    
    # ax1.axvline(255/30*0.93/127.5, lw = lw1*0.5, ls = 'dotted', dashes = [2, 2], color = 'black', alpha = 1.0, zorder = 10, label = r'$k_{D/L=0.118}$') # dashesで破線の間隔などを設定できる
    # ax1.axvline(255/40*0.93/127.5, lw = lw1*0.5, ls = 'dotted', color = 'black', alpha = 0.8, zorder = 10, label = r'$k_{D/L=0.157}$') # dashesで破線の間隔などを設定できる
    # ax1.axvline(255/50*0.93/127.5, lw = lw1*0.5, ls = 'dotted', dashes = [2, 2], color = 'magenta', alpha = 1.0, zorder = 10, label = r'$k_{D/L=0.196}$') # dashesで破線の間隔などを設定できる
    # ax1.axvline(255/60*0.93/127.5, lw = lw1*0.5, ls = 'dotted', dashes = [2, 2], color = 'magenta', alpha = 1.0, zorder = 10, label = r'$k_{D/L=0.235}$') # dashesで破線の間隔などを設定できる
    # ax1.axvline(255/70*0.93/127.5, lw = lw1*0.5, ls = 'dotted', dashes = [2, 2], color = 'deepskyblue', alpha = 1.0, zorder = 10, label = r'$k_{D/L=0.275}$') # dashesで破線の間隔などを設定できる
    
    # ax1.axvline(255/40*0.93/127.5, lw = lw1*0.5, ls = 'dotted', dashes = [1, 2], color = 'royalblue', label = r'$k_{D=40}$') # dashesで破線の間隔などを設定できる
    # za=np.arange(0.02, 1, 0.0001)
    # ua=0.7*10.0**(-6)*za**(-5/3)
    # ax1.plot(za, ua, lw = lw1*0.8, ls = 'dotted',  color = 'blue', 
    #         alpha = 1,     # 透明度
    #         clip_on = True, # プロット枠外にもプロットする
    #         zorder = 9,      # zorderが大きいほど前面に表示される
    #         label = r'$\propto k^{-5/3}$') 
    
    # 凡例の設定
    h1, l1 = ax1.get_legend_handles_labels()
    ax1.legend(h1, l1, 
    bbox_to_anchor = (1.0, 1.0), loc = "upper left", # bbox_to_anchorは凡例のlocの座標
    framealpha = 1.0, fancybox=False, fontsize=8.0,
    edgecolor = "black").get_frame().set_linewidth(alw*0.8)
    # ### 保存
    fig.savefig(plotdir + "we_d_break" + ext, bbox_inches = "tight") # bbox_inches="tight"で余白をなくす

##=================== main ===================##
if __name__ == '__main__':
    print("start main")
    print("datadir: ", datadir, ", plotdir: ", plotdir, ", ctxt: ", ctxt1)
    sns_set(fs1, tck_s1, alw, ctxt1)
    plot_y1y2()
    print("end main")