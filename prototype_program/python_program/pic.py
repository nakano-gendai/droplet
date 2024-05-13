import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys

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

fsize = 0   # 0：小さいサイズ（論文用）　1：大きいサイズ（プレゼン用）
ext = 'png' # 保存ファイルの拡張子　pdf,svg,pngなど
plotdir = 'plot_movie' # 保存用ディレクトリ
os.makedirs(plotdir, exist_ok = True) # ディレクトリ作成

datadir = '/home/nakano/anime/crit/re150/'
plotdir, datadir, ext = plotdir + '/', datadir + '/', '.' + ext

### サイズ、ラベルなどの設定
lx, ly = r'$x$', r'$z$' # r''でTeX文字にできる
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

if fsize == 0:
    fig = plt.figure(figsize = (4.0, 4.0), dpi = 100, linewidth = 0)
else:
    fig = plt.figure(figsize = (4.0, 4.0), dpi = 100, linewidth = 0)
ax1 = fig.add_subplot(111)
ax1.spines["top"].set_linewidth(alw)
ax1.spines["left"].set_linewidth(alw)
ax1.spines["bottom"].set_linewidth(alw)
ax1.spines["right"].set_linewidth(alw)

ax1.set_xlabel(lx, labelpad = lpad[0]) # 軸ラベル
ax1.set_ylabel(ly, labelpad = lpad[1])
xm, ym = [0, 299], [0, 101]
ax1.set_xlim(xm[0], xm[1]) # 軸の範囲
ax1.set_ylim(ym[0], ym[1])

print("start main")
print("datadir: ", datadir, ", plotdir: ", plotdir, ", ctxt: ", ctxt1)
sns_set(fs1, tck_s1, alw, ctxt1)

for n in range(33600,200000,100):
    if n >= 100 and n < 1000:
        dfile   = datadir + '0_{:03d}.d'.format(n)
    elif n >= 1000 and n < 10000:
        dfile   = datadir + '0_{:04d}.d'.format(n)
    elif n >= 10000 and n < 100000:
        dfile   = datadir + '0_{:05d}.d'.format(n)
    elif n >= 100000 and n < 1000000:
        dfile   = datadir + '0_{:06d}.d'.format(n)

    if os.path.exists(plotdir) == False or os.path.exists(datadir) == False: # plotdirやdatadirが存在しないときに中止する
        sys.exit("Error")

    ### 読み込み
    x1 = np.loadtxt(dfile, usecols = 2, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数int64など

    # uw=0.13
    # h=101
    # t = (n-5000)*2*uw/h
    # plt.title(r'$t^*={:2}$'.format(t), fontsize=10)
    a = x1.reshape(102,300)
    ax1.set_aspect("equal")
    extent1=[-0.5+xm[0], 0.5+xm[1], -0.5+ym[0], 0.5+ym[1]]
    im = ax1.imshow(a, origin="lower", cmap=plt.cm.plasma, vmin=2, vmax=5, extent=extent1)
    axpos = ax1.get_position()
    pp_ax = fig.add_axes([axpos.x1 + 0.02, axpos.y0, 0.05, axpos.height]) # 左下のx,y,幅,高さ?
    pp=fig.colorbar(im, ax=ax1, orientation="vertical", cax=pp_ax)
    pp.set_label(r"$\phi$", labelpad=5)


    ### 保存
    fig.savefig(plotdir + '{:08d}'.format(n) + ext, bbox_inches = "tight")

print("end main")



