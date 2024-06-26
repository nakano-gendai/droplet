#-------------------------------------------------------------------------------
# gnuplotの設定
#-------------------------------------------------------------------------------
reset
set nokey                 # 凡例の非表示
set xrange [0:80]       # x軸方向の範囲の設定
set yrange [0:80]       # y軸方向の範囲の設定
set size ratio -1           # 図を正方形にする
set pm3d map

set term gif animate      # 出力をgifアニメに設定
set output "taylortest2.gif"  # 出力ファイル名の設定

#-------------------------------------------------------------------------------
# 変数の設定
#-------------------------------------------------------------------------------
n0 = 1000    # ループ変数の初期値
n1 = 100000   # ループ変数の最大値
dn = 1000    # ループ変数の増加間隔

#-------------------------------------------------------------------------------
# ループの開始
#-------------------------------------------------------------------------------
load "anime.plt" 