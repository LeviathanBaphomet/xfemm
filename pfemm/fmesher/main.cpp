#include <iostream>
#include <string.h>
//#include "stdio.h"
#include "fmesher.h"

using std::string;

int main(int argc, char ** argv)
{

    FMesher MeshObj;
//    bool MeshUpToDate;
    string FilePath;
    //char tempFilePath[1024];

    if (argc < 2)
    {
        // request the file name from the user
        cout << "Enter fem file name:" << endl;
        getline(cin,FilePath);

        //scanf("%s", tempFilePath);
        // gets(tempFilePath);
        // FilePath = tempFilePath;
    }
    else if(argc > 2)
    {
        cout << "Too many input arguments" << endl;
        return -4;
    }
    else
    {
        FilePath = argv[1];
    }

    if (MeshObj.LoadFEMFile(FilePath) == false)
    {
        return -1;
    }

    if (MeshObj.HasPeriodicBC() == true)
    {
        if (MeshObj.DoPeriodicBCTriangulation(FilePath) != 0)
        {
            return -2;
        }
    }
    else
    {
        if (MeshObj.DoNonPeriodicBCTriangulation(FilePath) != 0)
        {
            return -3;
        }
    }

    //bool LoadMesh = MeshObj.LoadMesh(FilePath);

    //EndWaitCursor();

    //if(LoadMesh == true){

    //    //MeshUpToDate = true;

    //    //if(MeshFlag == false)
    //    //{
    //    //    OnShowMesh();
    //    //}
    //    //else
    //    //{
    //    //    InvalidateRect(NULL);
    //    //}

    //    //CStdString s;

    //    //s.Format("Created mesh with %i nodes", MeshObj.meshnode.size());

    //    //if (TheDoc->greymeshline.size() != 0)
    //    //    s += "\nGrey mesh lines denote regions\nthat have no block label.";
    //
    //    //if(bLinehook==false)
    //    //{
    //    //    WarnMessage(s,MB_ICONINFORMATION);
    //    //}
    //    //else
    //    //{
    //    //    lua_pushnumber(lua,TheDoc->meshnode.size());
    //    //}

    //}

    cout << "No errors" << endl;
    return 0;

}
