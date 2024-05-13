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
plotdir = 'plot3' # 保存用ディレクトリ
os.makedirs(plotdir, exist_ok = True) # ディレクトリ作成

# datadir = '/Users/nakanogendai/droplet/'
# datadir ='/home/nakano/anime/b4/data/data_fig6/'
datadir ="./"
plotdir, datadir, ext = plotdir + '/', datadir + '/', '.' + ext
dfile1   = datadir + 'taylor_para.d'
# dfile2   = datadir + 'ens_step_LES_medium.d'
# dfile3   = datadir + 'ens_step_LES_small.d'
# dfile4   = datadir + 'test64.d'
# dfile5   = datadir + 'ens_step_LES_medium.d'
# dfile6   = datadir + 'ens_step_LES_small.d'
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
x1 = np.loadtxt(dfile1, usecols = 1, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
# x2 = np.loadtxt(dfile2, usecols = 1, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
# x3 = np.loadtxt(dfile3, usecols = 1, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
# x4 = np.loadtxt(dfile4, usecols = 0, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
# x5 = np.loadtxt(dfile5, usecols = 1, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
# x6 = np.loadtxt(dfile6, usecols = 1, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
# x7 = np.loadtxt(dfile7, usecols = 0, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数
# x8 = np.loadtxt(dfile8, usecols = 0, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数

y1 = np.loadtxt(dfile1, usecols = 3, dtype = 'float64')
# y2 = np.loadtxt(dfile2, usecols = 2, dtype = 'float64')
# y3 = np.loadtxt(dfile3, usecols = 2, dtype = 'float64')
# y4 = np.loadtxt(dfile4, usecols = 1, dtype = 'float64')
# y5 = np.loadtxt(dfile5, usecols = 2, dtype = 'float64')
# y6 = np.loadtxt(dfile6, usecols = 2, dtype = 'float64')
# y7 = np.loadtxt(dfile7, usecols = 3, dtype = 'float64')
# y8 = np.loadtxt(dfile8, usecols = 3, dtype = 'float64')



### サイズ、ラベルなどの設定
lx, ly = r'$tU_{\mathrm{max}}/D_{\mathrm{v}}$', r'$\mathrm{Re}_\lambda$' # r''でTeX文字にできる
# lx, ly = r'$tU_{\mathrm{max}}/D_{\mathrm{v}}$', r'$\mathrm{Re}_{\lambda}$' # r''でTeX文字にできる
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
    # xm, ym = [1500, 1950], [0.00000000001, 0.00000000002]
    xm, ym = [1500, 1950], [110,150]
    ax1.set_xlim(xm[0], xm[1]) # 軸の範囲
    ax1.set_ylim(ym[0], ym[1])
    # ax1.set_yticks(np.arange(0.0, ym[1] + 0.0000005)) # xmaxまで0.2刻みの目盛り線
    ax1.tick_params(axis='x', pad = tpad[0])
    ax1.tick_params(axis='y', pad = tpad[1])
    # plt.yscale('log')
    # plt.xscale('log')

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
    # ax1.plot(x1, y1, lw = 1, ls = 'none', marker = 'o', ms = ms1*0.4, color = 'black', 
    #         alpha = 1,     # 透明度
    #         clip_on = False, # プロット枠外にもプロットする
    #         zorder = 2,
    #         label = 'break') 
    # ax1.plot(x1, y1, lw = 1, ls = 'none', marker = 'o', ms = ms1*2.5, mew = lw1*1.5, mfc = 'none', color = 'blue', 
    #         alpha = 1,     # 透明度
    #         clip_on = False, # プロット枠外にもプロットする
    #         zorder = 1,      # zorderが大きいほど前面に表示される
    #         label = 'probability density function') 
    # ax1.plot(x3, y3, lw = 1, ls = 'solid', marker = '^', ms = ms1*1.0, mew = lw1*1.0, mfc = 'none', color = 'tomato', 
    #         alpha = 1,     # 透明度
    #         clip_on = False, # プロット枠外にもプロットする
    #         zorder = 1,      # zorderが大きいほど前面に表示される
    #         label = r'$\mathrm{Re}\approx 0.2, \mathrm{Ca}\approx 0.3$') 
    # ax1.plot(x4, y4, lw = 1, ls = 'solid', marker = 's', ms = ms1*1.0, mew = lw1*1.0, mfc = 'none', color = 'royalblue', 
    #         alpha = 1,     # 透明度
    #         clip_on = False, # プロット枠外にもプロットする
    #         zorder = 1,      # zorderが大きいほど前面に表示される
    #         label = r'$\mathrm{Re}\approx 0.2, \mathrm{Ca}\approx 0.4$') 
    # ax1.plot(x5, y5, lw = 1, ls = 'solid', marker = '+', ms = ms1*1.0, mew = lw1*1.0, mfc = 'none', color = 'lime', 
    #         alpha = 1,     # 透明度
    #         clip_on = False, # プロット枠外にもプロットする
    #         zorder = 2,      # zorderが大きいほど前面に表示される
    #         label = r'$\mathrm{Re}\approx 0.2$, $\mathrm{Ca}\approx 0.1$ (initial data)') 
    # ax1.plot(x6, y6, lw = 1, ls = 'solid', marker = '+', ms = ms1*1.0, mew = lw1*1.0, mfc = 'none', color = 'tomato', 
    #         alpha = 1,     # 透明度
    #         clip_on = False, # プロット枠外にもプロットする
    #         zorder = 2,      # zorderが大きいほど前面に表示される
    #         label = r'$\mathrm{Re}\approx 0.2$, $\mathrm{Ca}\approx 0.2$ (initial data)') 
    # ax1.plot(x7, y7, lw = 0.1, ls = 'solid', marker = 'D', ms = ms1*1.0, mew = lw1*1.0, mfc = 'none', color = 'lime', 
    #         alpha = 1,     # 透明度
    #         clip_on = False, # プロット枠外にもプロットする
    #         zorder = 2,      # zorderが大きいほど前面に表示される
    #         label = r'$\mathrm{Re}\approx 0.2$, $\mathrm{Ca}\approx 0.3$ (initial data)') 
    # ax1.plot(x8, y8, lw = 0.1, ls = 'solid', marker = '*', ms = ms1*1.0, mew = lw1*1.0, mfc = 'none', color = 'gray', 
    #         alpha = 1,     # 透明度
    #         clip_on = False, # プロット枠外にもプロットする
    #         zorder = 2,      # zorderが大きいほど前面に表示される
    #         label = r'$\mathrm{Re}\approx 0.2$, $\mathrm{Ca}\approx 0.4$ (initial data)') 

    
    ax1.plot(x1, y1, lw = lw1*1.5, ls = 'solid',  color = 'black', alpha = 1.0, clip_on = True, zorder = 13)
    # ax1.plot(x2, y2, lw = lw1*1.5, ls = 'dotted',  color = 'blue', alpha = 1.0, clip_on = True, zorder = 15, label = r"medium (LES, $\Delta=8\Delta x$)")
    # ax1.plot(x3, y3, lw = lw1*1.5, ls = 'dashed',  color = 'orange', alpha = 1.0, clip_on = True, zorder = 16, label = r"small (LES, $\Delta=8\Delta x$)")
    # ax1.plot(x4, y4, lw = lw1*1.2, ls = 'dashdot',  color = 'orange', alpha = 1.0, clip_on = True, zorder = 17, label = r'LES ($64^3$)')
    # ax1.plot(x5, y5, lw = lw1*1.5, ls = 'dotted',  color = 'blue', alpha = 1, clip_on = True, zorder = 11, label = 'medium (LES)')
    # ax1.plot(x6, y6, lw = lw1*1.5, ls = 'dashed',  color = 'blue', alpha = 1, clip_on = True, zorder = 11, label = 'small (LES)')
    # ax1.plot(x7, y7, lw = lw1*2.3, ls = 'dashed',  color = 'cyan', alpha = 1, clip_on = True, zorder = 11, label = r'$\mathrm{Re}\approx 0.2$, $\mathrm{Ca}\approx 0.3$ (initial data)')
    # ax1.plot(x8, y8, lw = lw1*2.3, ls = 'solid',  color = 'cyan', alpha = 1, clip_on = True, zorder = 11, label = r'$\mathrm{Re}\approx 0.2$, $\mathrm{Ca}\approx 0.4$ (initial data)')
    
    # ax1.axvline(0.392, lw = lw1*0.8, ls = 'dotted', dashes = [1, 2], color = 'red', label = r'$\pi / \Delta$') # dashesで破線の間隔などを設定できる
    # za=np.arange(0.02, 1, 0.0001)
    # ua=0.7*10.0**(-6)*za**(-5/3)
    # ax1.plot(za, ua, lw = lw1*0.8, ls = 'dotted',  color = 'blue', 
    #         alpha = 1,     # 透明度
    #         clip_on = True, # プロット枠外にもプロットする
    #         zorder = 9,      # zorderが大きいほど前面に表示される
    #         label = r'$\propto k^{-5/3}$') 
    
    # 凡例の設定
    # h1, l1 = ax1.get_legend_handles_labels()
    # ax1.legend(h1, l1, 
    # bbox_to_anchor = (1.0, 1.0), loc = "upper left", # bbox_to_anchorは凡例のlocの座標
    # framealpha = 1.0, fancybox=False, fontsize=8.0,
    # edgecolor = "black").get_frame().set_linewidth(alw*0.8)
    ### 保存
    fig.savefig(plotdir + "taylor_lambda_256_re80000" + ext, bbox_inches = "tight") # bbox_inches="tight"で余白をなくす

##=================== main ===================##
if __name__ == '__main__':
    print("start main")
    print("datadir: ", datadir, ", plotdir: ", plotdir, ", ctxt: ", ctxt1)
    sns_set(fs1, tck_s1, alw, ctxt1)
    plot_y1y2()
    print("end main")