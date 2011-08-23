function [ output ] = DVImgGetExtInfo( Stream )

if Stream == 0
    helpdlg('Retrieves information about the metadata fields in the image''s extended header. A struct is returned, with fields corresponding to each entry in the extended header and the value of each field being the data type of that field. See DVImgGetExtHeaderField and DVImgSetExtHeaderField.','int DVImgGetExtInfo(int Stream)');
else
    NumFields = calllib(DVImgLibName,'DVImgGetNumFields',Stream);
    for Field=0:NumFields-1
        FieldName=calllib(DVImgLibName,'DVImgGetFieldName',Stream,Field);
        FieldType=char(calllib(DVImgLibName,'DVImgGetFieldType',Stream,Field));
        Fields.(FieldName)={FieldType};
    end

    if NumFields == 0
        Fields = struct();
    end

    output = Fields;
end 