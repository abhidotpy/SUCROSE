import os
import subprocess

p1 = subprocess.run(["ifx", "-r8", "-Ofast", "-xHost", "-qipo", "-qopenmp", "-module", "modules", "-c", "fortran/xmath.f90", "-o", "build/xmath.o"], capture_output=True)

if (p1.returncode == 0):
    print("xmath.f90 build successful")
else:
    print("xmath.f90 build failed to run")
    print(p1.stderr)
    exit()

p2 = subprocess.run(["ifx", "-r8", "-Ofast", "-xHost", "-qipo", "-qopenmp", "-module", "modules", "-c", "fortran/system.f90", "-o", "build/system.o"], capture_output=True)

if (p2.returncode == 0):
    print("system.f90 build successful")
else:
    print("system.f90 build failed to run")
    print(p2.stderr)
    exit()

p3 = subprocess.run(["ifx", "-r8", "-Ofast", "-xHost", "-qipo", "-qopenmp", "-module", "modules", "-c", "fortran/force.f90", "-o", "build/force.o"], capture_output=True)

if (p3.returncode == 0):
    print("force.f90 build successful")
else:
    print("force.f90 build failed to run")
    print(p3.stderr)
    exit()

p4 = subprocess.run(["ifx", "-r8", "-Ofast", "-xHost", "-qipo", "-qopenmp", "-module", "modules", "-c", "fortran/integrators.f90", "-o", "build/integrators.o"], capture_output=True)

if (p4.returncode == 0):
    print("integrators.f90 build successful")
else:
    print("integrators.f90 build failed to run")
    print(p4.stderr)
    exit()

p5 = subprocess.run(["ifx", "-r8", "-Ofast", "-xHost", "-qipo", "-qopenmp", "-module", "modules", "-c", "fortran/reporters.f90", "-o", "build/reporters.o"], capture_output=True)

if (p5.returncode == 0):
    print("reporters.f90 build successful")
else:
    print("reporters.f90 build failed to run")
    print(p5.stderr)
    exit()

