from subprocess import run, Popen, PIPE, CREATE_NEW_CONSOLE
import os

#subprocess.call(env, cwd=prjDir, shell=True)

prjDir = "d:/proj/diffusion_solver/"
my_env = os.environ.copy()
my_env["PATH"] = "c:/proj/adolc/ADOL-C-2.6.3/ADOL-C-2.6.3/MSVisualStudio/v14/ColPack/Release/;C:/Program Files/VTK/bin;c:/proj/adolc/ADOL-C-2.6.3/ADOL-C-2.6.3/MSVisualStudio/v14/x64/sparse/;" + my_env["PATH"] 
cmd = 'x64\Release\diffusion_solver.exe --dir snaps/p250_k500_av_len/ --pwf 250 --perm 500.0 --calc_duration 0.25 --init_step 0.000005 --frac_len 120.0 \n'
p = Popen(cmd, cwd=prjDir, env=my_env, encoding='cp866')
