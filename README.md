# Spectra for Visual Studio

This is a Visual Studio solution for [Spectra](https://github.com/yixuan/spectra).

## Instructions

The repository contains a C API for the Spectra eigensolver. The project depends on the currently (2020-06-30) unreleased 1.y.z branch of [Spectra](https://github.com/yixuan/spectra/archive/1.y.z.zip). Spectra depends on [Eigen](https://gitlab.com/libeigen/eigen). Make sure that both libraries are available in the include subfolder.

Pre-compiled binaries for windows users can be found [here](http://wo80.bplaced.net/math/packages.html).

## Why?

The project was created to maintain ARPACK builds matching the [CSparse.Interop](https://github.com/wo80/csparse-interop) bindings for C#.