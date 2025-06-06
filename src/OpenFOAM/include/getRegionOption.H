/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Get the value from the -region option

Required Classes
    - Foam::polyMesh

Required Variables
    - args [argList]

Provided Variables
    - regionName [word]

\*---------------------------------------------------------------------------*/

Foam::word regionName
(
    args.getOrDefault<Foam::word>("region", Foam::polyMesh::defaultRegion)
);


// ************************************************************************* //
