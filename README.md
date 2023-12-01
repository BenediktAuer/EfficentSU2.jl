# EfficentSU2

[![Build Status](https://github.com/BenediktAuer/EfficentSU2.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/BenediktAuer/EfficentSU2.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![CompatHelper](https://github.com/BenediktAuer/EfficentSU2.jl/actions/workflows/CompatHelper.yml/badge.svg)](https://github.com/BenediktAuer/EfficentSU2.jl/actions/workflows/CompatHelper.yml)

For iterating over even lattice sites Iterators.filter(x->(-1)^sum(x.I)>0,indices) where indices = CartesianIndices(matrix)
for the odd use Iterators.filter(x->(-1)^sum(x.I)<0,indices) these can be paralizied because every odd site only even sites are needed
