const polyMesh& mesh = this->internalField().mesh()();

const volScalarField& field =
            this->db().objectRegistry::lookupObject<volScalarField>("C");
const volScalarField& auxField =
            this->db().objectRegistry::lookupObject<volScalarField>("Deff");
const label& patchID = this->patch().index();
const IOdictionary& IOd
          = this->db().lookupObject<IOdictionary>("transportProperties");
scalar scalarVal =  (new dimensionedScalar("lR", dimLength, IOd))->value();
scalar scalarAux =  (new dimensionedScalar("D", (dimArea/dimTime), IOd))->value();
scalarField gradField = -field.boundaryField()[patchID].snGrad();
vectorField faceNorm = mesh.boundaryMesh()[patchID].faceNormals();
    
const scalar dt = this->db().time().deltaTValue();
vectorField fD = dt * scalarVal * gradField * faceNorm;
fD = fD * auxField.boundaryField()[patchID] / scalarAux;

this->setDisp(fD);

Info << "Tohoko - libs/boundaryConditions/gradCMotionDeff.H" << endl;
