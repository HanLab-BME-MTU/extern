function [ output ] = DVImgGetExtHeaderField( Stream,FieldName )

if Stream == 0
    helpdlg('Retrieves a metadata field from the extended header. If the field does not exist an error is thrown. See DVImgSetExtHeaderField and DVImgGetExtInfo.','int DVImgGetHeaderField(int Stream, char FieldName)');
else
    Fields = DVImgGetExtInfo(Stream);
    FieldType = char(Fields.(FieldName)(1));

    switch FieldType
        case 'd'
            output = calllib(DVImgLibName,'DVImgGetFieldDouble',Stream,FieldName);
        case 'f'
            output = calllib(DVImgLibName,'DVImgGetFieldFloat',Stream,FieldName);
        case 'i'
            output = calllib(DVImgLibName,'DVImgGetFieldInt',Stream,FieldName);
        case 'c'
            output = calllib(DVImgLibName,'DVImgGetFieldChar',Stream,FieldName);
        otherwise
            error('Field type must be ''d'', ''f'', ''i'', or ''c''');
    end

end
