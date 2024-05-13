"""
Line plot etc. LaTeX font
10/01/2023
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
# import matplotlib.gridspec as gridspec

### 入出力の設定
fsize = 0   # 0：小さいサイズ（論文用）　1：大きいサイズ（プレゼン用）
ext = 'pdf' # 保存ファイルの拡張子　pdf,svg,pngなど
plotdir = 'plot' # 保存用ディレクトリ
os.makedirs(plotdir, exist_ok = True) # ディレクトリ作成

datadir = '.'
# datadir = '/home/nakano/droplet/program/eta0.1re10ca0.5break/'
plotdir, datadir, ext = plotdir + '/', datadir + '/', '.' + ext
dfile1   = datadir + '5800.d'
dfile2   = datadir + '18000.d'
dfile3   = datadir + '21200.d'
dfile4   = datadir + '9400.d'
dfile5   = datadir + '12100.d'
dfile6   = datadir + '13100.d'

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

#figure()でグラフを表示する領域をつくり，figというオブジェクトにする．
fig = plt.figure()
plt.subplots_adjust(wspace=0.08, hspace=-0.4)

# add_subplot()でグラフを描画する領域を追加する．引数は行，列，場所
# fig = plt.figure(figsize = (3.0, 3.0))
ax1 = fig.add_subplot(3, 2, 1)
ax2 = fig.add_subplot(3, 2, 3)
ax3 = fig.add_subplot(3, 2, 5)
ax4 = fig.add_subplot(3, 2, 2)
ax5 = fig.add_subplot(3, 2, 4)
ax6 = fig.add_subplot(3, 2, 6)



### 読み込み
x1 = np.loadtxt(dfile1, usecols = 6, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数int64など
x2 = np.loadtxt(dfile2, usecols = 6, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数int64など
x3 = np.loadtxt(dfile3, usecols = 6, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数int64など
x4 = np.loadtxt(dfile4, usecols = 7, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数int64など
x5 = np.loadtxt(dfile5, usecols = 7, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数int64など
x6 = np.loadtxt(dfile6, usecols = 7, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数int64など

# y1 = np.loadtxt(dfile, usecols = 0, dtype = 'float64')
# y2 = np.loadtxt(dfile, usecols = 3, dtype = 'float64')
# y3 = np.loadtxt(dfile, usecols = 3, dtype = 'float64')

### サイズ、ラベルなどの設定
lx, ly = r"$x / H$", r"$z / H$" # r''でTeX文字にできる
if fsize == 0:
    fs1, lw1, ms1 = 1., 1., 2.
    tck_s1, alw = 3, 0.625
    ctxt1 = 'paper'
    lpad = [5, 5] # 軸とラベルの間隔
    tpad = [5, 5] # 軸と数値の間隔
else:
    fs1, lw1, ms1 = 1.5, 3., 14.
    tck_s1, alw = 7, 1.25
    ctxt1, ext = 'talk', '_talk' + ext
    lpad = [10, 8]
    tpad = [8, 12]


def plot_y1y2(): # y1, y2プロット用
    # if fsize == 0:
    #     fig = plt.figure(figsize = (2.3, 2.3), dpi = 100, linewidth = 0)
    # else:
    #     fig = plt.figure(figsize = (7.2, 4.8), dpi = 100, linewidth = 0)
    # ax1 = fig.add_subplot(111)
    
    # ax1.spines["top"].set_linewidth(alw)
    # ax1.spines["left"].set_linewidth(alw)
    # ax1.spines["bottom"].set_linewidth(alw)
    # ax1.spines["right"].set_linewidth(alw)
    
    # ax2.spines["top"].set_linewidth(alw)
    # ax2.spines["left"].set_linewidth(alw)
    # ax2.spines["bottom"].set_linewidth(alw)
    # ax2.spines["right"].set_linewidth(alw)
    
    # ax3.spines["top"].set_linewidth(alw)
    # ax3.spines["left"].set_linewidth(alw)
    # ax3.spines["bottom"].set_linewidth(alw)
    # ax3.spines["right"].set_linewidth(alw)
    
    # ax4.spines["top"].set_linewidth(alw)
    # ax4.spines["left"].set_linewidth(alw)
    # ax4.spines["bottom"].set_linewidth(alw)
    # ax4.spines["right"].set_linewidth(alw)
    
    # ax5.spines["top"].set_linewidth(alw)
    # ax5.spines["left"].set_linewidth(alw)
    # ax5.spines["bottom"].set_linewidth(alw)
    # ax5.spines["right"].set_linewidth(alw)
    
    # ax6.spines["top"].set_linewidth(alw)
    # ax6.spines["left"].set_linewidth(alw)
    # ax6.spines["bottom"].set_linewidth(alw)
    # ax6.spines["right"].set_linewidth(alw)

    ax1.tick_params(labelbottom=False, labelleft=True, labelright=False, labeltop=False)
    ax2.tick_params(labelbottom=False, labelleft=True, labelright=False, labeltop=False)
    ax3.tick_params(labelbottom=True, labelleft=True, labelright=False, labeltop=False)
    ax4.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False)
    ax5.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False)
    ax6.tick_params(labelbottom=True, labelleft=False, labelright=False, labeltop=False)

    # plt.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False, bottom=False, left=False, right=False, top=False)

    # ax1.set_xlabel(lx, labelpad = lpad[0]) # 軸ラベル
    ax1.set_ylabel(ly, labelpad = lpad[1])
    # ax2.set_xlabel(lx, labelpad = lpad[0]) # 軸ラベル
    ax2.set_ylabel(ly, labelpad = lpad[1])
    ax3.set_xlabel(lx, labelpad = lpad[0]) # 軸ラベル
    ax3.set_ylabel(ly, labelpad = lpad[1])
    # ax4.set_xlabel(lx, labelpad = lpad[0]) # 軸ラベル
    # ax4.set_ylabel(ly, labelpad = lpad[1])
    # ax5.set_xlabel(lx, labelpad = lpad[0]) # 軸ラベル
    # ax5.set_ylabel(ly, labelpad = lpad[1])
    ax6.set_xlabel(lx, labelpad = lpad[0]) # 軸ラベル
    # ax6.set_ylabel(ly, labelpad = lpad[1])
    

    ax1.tick_params(axis='x', pad = tpad[0])
    ax1.tick_params(axis='y', pad = tpad[1])
    xm1, ym1 = [0, 299/71], [0, 71/71]
    ax1.set_xlim(xm1[0], xm1[1]) # 軸の範囲
    ax1.set_ylim(ym1[0], ym1[1])
    ax2.set_xlim(xm1[0], xm1[1]) # 軸の範囲
    ax2.set_ylim(ym1[0], ym1[1])
    ax3.set_xlim(xm1[0], xm1[1]) # 軸の範囲
    ax3.set_ylim(ym1[0], ym1[1])
    xm2, ym2 = [0, 299/71], [0, 71/71]
    ax4.set_xlim(xm2[0], xm2[1]) # 軸の範囲
    ax4.set_ylim(ym2[0], ym2[1])
    ax5.set_xlim(xm2[0], xm2[1]) # 軸の範囲
    ax5.set_ylim(ym2[0], ym2[1])
    ax6.set_xlim(xm2[0], xm2[1]) # 軸の範囲
    ax6.set_ylim(ym2[0], ym2[1])

    # ax1.set_xticks(0,299/71,1)
    # ax1.set_yticks(0,71/71,0.5)
    
    # ax1.set_xticks(np.arange(0, xm[1] + 1.e-3, )) # xmaxまで0.2刻みの目盛り線
    # ax1.tick_params(axis='x', pad = tpad[0])
    # ax1.tick_params(axis='y', pad = tpad[1])

    ### 水平線と鉛直線
    # ax1.axhline(0, lw = lw1*0.8, ls = 'dotted', label = r'$y=0$')
    # ax1.axvline(4.8, lw = lw1*0.8, ls = 'dashed', dashes = [2, 5], color = 'C2', label = r'$x=x\ssmr{c}$') # dashesで破線の間隔などを設定できる
    # ytext = -22.
    # ax1.axhline(ytext, lw = lw1*0.8, c = 'C4', ls = 'dashdot')
    # ax1.text(0.7, (ytext - ym[0])/(ym[1] - ym[0]) + 0.05, r'$y=%4.1f$'%(ytext), transform = ax1.transAxes) # 図中にtextを入れる

    # za=np.arange(-1,1.1,0.001)
    # ua=2*za-1
    # ax1.plot(za, ua, lw = lw1*0.5, ls = 'solid',  color = 'black', 
    #         alpha = 1,     # 透明度
    #         clip_on = True, # プロット枠外にもプロットする
    #         zorder = 9,      # zorderが大きいほど前面に表示される
    #         label = 'Theory') 
    
    bin=0.5/71

    a1 = x1.reshape(72,300)
    ax1.set_aspect("equal")
    extent1=[-bin+xm1[0], bin+xm1[1], -bin+ym1[0], bin+ym1[1]]
    im = ax1.imshow(a1, origin="lower", cmap=plt.cm.plasma, vmin=2, vmax=5, extent=extent1)
    # axpos = ax1.get_position()
    # pp_ax = fig.add_axes([axpos.x1 + 0.01, axpos.y0, 0.015, axpos.height]) # 左下のx,y,幅,高さ?
    # pp=fig.colorbar(im, ax=ax1, orientation="vertical", cax=pp_ax)
    # pp.set_label(r"$\phi$", labelpad=0.2)
    
    a2 = x2.reshape(72,300)
    ax2.set_aspect("equal")
    extent2=[-bin+xm1[0], bin+xm1[1], -bin+ym1[0], bin+ym1[1]]
    im = ax2.imshow(a2, origin="lower", cmap=plt.cm.plasma, vmin=2, vmax=5, extent=extent2)
    # axpos = ax2.get_position()
    # pp_ax = fig.add_axes([axpos.x1 + 0.01, axpos.y0, 0.015, axpos.height]) # 左下のx,y,幅,高さ?
    # pp=fig.colorbar(im, ax=ax2, orientation="vertical", cax=pp_ax)
    # pp.set_label(r"$\phi$", labelpad=0.2)
    
    a3 = x3.reshape(72,300)
    ax3.set_aspect("equal")
    extent3=[-bin+xm1[0], bin+xm1[1], -bin+ym1[0],bin+ym1[1]]
    im = ax3.imshow(a3, origin="lower", cmap=plt.cm.plasma, vmin=2, vmax=5, extent=extent3)
    # axpos = ax3.get_position()
    # pp_ax = fig.add_axes([axpos.x1 + 0.01, axpos.y0, 0.015, axpos.height]) # 左下のx,y,幅,高さ?
    # pp=fig.colorbar(im, ax=ax3, orientation="vertical", cax=pp_ax)
    # pp.set_label(r"$\phi$", labelpad=0.2)
    
    a4 = x4.reshape(66,255)
    ax4.set_aspect("equal")
    extent4=[-bin+xm2[0], bin+xm2[1], -bin+ym2[0], bin+ym2[1]]
    im = ax4.imshow(a4, origin="lower", cmap=plt.cm.plasma, vmin=2, vmax=5, extent=extent4)
    # axpos = ax4.get_position()
    # pp_ax = fig.add_axes([axpos.x1 + 0.01, axpos.y0, 0.015, axpos.height]) # 左下のx,y,幅,高さ?
    # pp=fig.colorbar(im, ax=ax4, orientation="vertical", cax=pp_ax)
    # pp.set_label(r"$\phi$", labelpad=0.2)
    
    a5 = x5.reshape(66,255)
    ax5.set_aspect("equal")
    extent5=[-bin+xm2[0], bin+xm2[1], -bin+ym2[0], bin+ym2[1]]
    im = ax5.imshow(a5, origin="lower", cmap=plt.cm.plasma, vmin=2, vmax=5, extent=extent5)
    # axpos = ax5.get_position()
    # pp_ax = fig.add_axes([axpos.x1 + 0.01, axpos.y0, 0.015, axpos.height]) # 左下のx,y,幅,高さ?
    # pp=fig.colorbar(im, ax=ax5, orientation="vertical", cax=pp_ax)
    # pp.set_label(r"$\phi$", labelpad=0.2)
    
    a6 = x6.reshape(66,255)
    ax6.set_aspect("equal")
    extent6=[-bin+xm2[0], bin+xm2[1], -bin+ym2[0], bin+ym2[1]]
    im = ax6.imshow(a6, origin="lower", cmap=plt.cm.plasma, vmin=2, vmax=5, extent=extent6)
    # axpos = ax6.get_position()
    # pp_ax = fig.add_axes([axpos.x1 + 0.01, axpos.y0, 0.015, axpos.height]) # 左下のx,y,幅,高さ?
    # pp=fig.colorbar(im, ax=ax6, orientation="vertical", cax=pp_ax)
    # pp.set_label(r"$\phi$", labelpad=0.2)

    #colorbarを1つにまとめる
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.83, 0.25, 0.015, 0.5])
    pp = fig.colorbar(im, cax=cbar_ax)
    pp.set_label(r"$\phi$", labelpad=0.4)
    
    # fig.tight_layout()  
    # plt.colorbar()
    # ### 読み込みデータのプロット
    # ax1.plot(x1, y1, lw = 0, ls = 'solid', marker = 'o', ms = ms1*1.7, mew = lw1*0.5, mfc = 'none', color = 'C0', 
    #         alpha = 1,     # 透明度
    #         clip_on = False, # プロット枠外にもプロットする
    #         zorder = 8,      # zorderが大きいほど前面に表示される
    #         label = 'LBM') 
    # ax1.plot(x1, y2, lw = lw1*0.5, ls = 'solid', color = 'black', alpha = 1, clip_on = True, zorder = 9, label = 'Theory')

    # ### 凡例の設定
    # h1, l1 = ax1.get_legend_handles_labels()
    # ax1.legend(h1, l1, 
    # bbox_to_anchor = (1.0, 1.0), loc = "upper left", # bbox_to_anchorは凡例のlocの座標
    # framealpha = 1.0, fancybox=False, 
    # edgecolor = "black").get_frame().set_linewidth(alw*0.8)
    ### 保存
    fig.savefig(plotdir + "test2" + ext, bbox_inches = "tight") # bbox_inches="tight"で余白をなくす

##=================== main ===================##
if __name__ == '__main__':
    print("start main")
    print("datadir: ", datadir, ", plotdir: ", plotdir, ", ctxt: ", ctxt1)
    sns_set(fs1, tck_s1, alw, ctxt1)
    plot_y1y2()
    print("end main")