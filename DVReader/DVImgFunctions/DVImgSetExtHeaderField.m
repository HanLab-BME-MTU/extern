function [ output ] = DVImgSetExtHeaderField( Stream,FieldName,FieldType,FieldValue )

if Stream == 0
    helpdlg('Writes to a metadata field in the extended header. If the field already exists it will be overwritten, if not, it will be created. FieldType is a single char indicating the datatype that the field will be written in -- ''c'' for char, ''f'' for float, ''i'' for int, or ''d'' for double.  See DVImgGetExtHeaderField and DVImgGetExtInfo.','int DVImgSetHeaderField(int Stream, char FieldName, char FieldType, FieldValue)');
else
    switch FieldType
        case 'd'
            output = calllib(DVImgLibName,'DVImgSetFieldDouble',Stream,FieldName,FieldValue);
        case 'f'
            output = calllib(DVImgLibName,'DVImgSetFieldFloat',Stream,FieldName,FieldValue);
        case 'i'
            output = calllib(DVImgLibName,'DVImgSetFieldInt',Stream,FieldName,FieldValue);
        case 'c'
            output = calllib(DVImgLibName,'DVImgSetFieldChar',Stream,FieldName,FieldValue);
        otherwise
            error('Field type must be ''d'', ''f'', ''i'', or ''c''');
    end

    DVImgPrintErrText(output);
end