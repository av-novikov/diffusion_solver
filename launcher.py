from subprocess import run, Popen, PIPE, CREATE_NEW_CONSOLE
import os

#subprocess.call(env, cwd=prjDir, shell=True)

prjDir = "f:/Temperature/diffusion_solver/"
my_env = os.environ.copy()
my_env["PATH"] = "C:/Program Files/VTK/bin;f:/Temperature/ADOL-C-2.6.2/ADOL-C-2.6.3/MSVisualStudio/v14/x64/sparse/;" + my_env["PATH"] 

cmd1 = 'x64\Release\diffusion_solver.exe --dir snaps/p240_k50_av_len_120/ --pwf 240 --perm 50.0 --calc_duration 1.0 --init_step 0.000005 --frac_len 120.0 \n'
cmd2 = 'x64\Release\diffusion_solver.exe --dir snaps/p250_k50_av_len_120/ --pwf 250 --perm 50.0 --calc_duration 1.0 --init_step 0.000005 --frac_len 120.0 \n'
cmd3 = 'x64\Release\diffusion_solver.exe --dir snaps/p260_k50_av_len_120/ --pwf 260 --perm 50.0 --calc_duration 0.7 --init_step 0.000005 --frac_len 120.0 \n'
cmd4 = 'x64\Release\diffusion_solver.exe --dir snaps/p280_k50_av_len_120/ --pwf 280 --perm 50.0 --calc_duration 0.5 --init_step 0.000005 --frac_len 120.0 \n'
cmd5 = 'x64\Release\diffusion_solver.exe --dir snaps/p280_k50_av_recalc/ --pwf 280 --perm 50.0 --calc_duration 0.5 --init_step 0.000005 --frac_len 240.0 \n'

#cmd_debug = 'x64/Debug/diffusion_solver.exe --dir snaps/ --pwf 260 --perm 50 \n'
p1 = Popen(cmd1, cwd=prjDir, env=my_env, encoding='cp866')
p2 = Popen(cmd2, cwd=prjDir, env=my_env, encoding='cp866')
p3 = Popen(cmd3, cwd=prjDir, env=my_env, encoding='cp866')
p4 = Popen(cmd4, cwd=prjDir, env=my_env, encoding='cp866')
p4 = Popen(cmd5, cwd=prjDir, env=my_env, encoding='cp866')
