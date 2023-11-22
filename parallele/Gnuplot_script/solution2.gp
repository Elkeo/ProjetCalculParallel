set palette defined (-5 0 0 1, 0 1 1 1, 5 1 0 0)
set terminal gif enhanced font Arial 30 animate delay 50 loop 1 optimize size 2000,2000
set output "solution2/solution2.gif"
set xr[0:1]
set yr[0:1]
dt=0.1
do for [i=0:11] {
t=i*dt
set title "t = ".sprintf("%f", t)." s"
plot "solution2/solutionFile_".i.".dat" u 1:2:3 with image
}
set output