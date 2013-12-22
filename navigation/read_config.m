function user = read_config( filename )

user = struct;
user.batch_date = datestr( now, 30 );

fid = fopen( filename );

curr_item = [];
cnt = 0;
while( 1 )
  tline = fgetl( fid );
  if( ~ischar( tline ) )
    break;
  end
  cnt = cnt + 1;

  tline = strtrim( tline );
  if( length( tline ) == 0 || tline(1) == '%' ) % line is empty or a comment
    continue;
  end

  if( tline(1) == '[' && tline(end) == ']' ) % line is an item name
    curr_item = strtrim( tline(2:end-1) );
    continue;
  end

  idx = find( tline == '=', 1 );
  if( ~isempty( idx ) )
    name = strtrim( tline( 1:idx-1 ) );
    value = strtrim( tline( idx+1:end) );
    if( strcmp( curr_item, 'main' ) == 1 )
      str = sprintf( 'user.%s = %s;', name, value );
      eval( str );
    else
      str = sprintf( 'user.%s.%s = %s;', curr_item, name, value );
      eval( str );
    end
    continue;
  end  

  error( 'Could not process line %d: %s', cnt, tline );
end

status = fclose( fid );
if( status ~= 0 )
  error( 'Closing file %s finished with an error.', filename );
end

%fprintf( 1, 'Successfully read %d lines from configuration file %s.\n', cnt, filename );
