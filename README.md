# Spectra for Visual Studio

This is a Visual Studio solution for [Spectra](https://github.com/yixuan/spectra).

## Instructions

The repository contains a C API for the Spectra eigensolver. The project depends on the currently (2021-02-17) unreleased master branch of [Spectra](https://github.com/yixuan/spectra/archive/master.zip). Spectra depends on [Eigen](https://gitlab.com/libeigen/eigen). Make sure that both libraries are available in the include subfolder.

Pre-compiled binaries for windows users can be found [here](http://wo80.bplaced.net/math/packages.html).

## Why?

The project was created to maintain Spectra builds matching the [CSparse.Interop](https://github.com/wo80/csparse-interop) bindings for C#.
