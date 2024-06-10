# This is the PNG module (part of MeiSai Framework)

# Loads Zlib library for (de-)compression of the image data
require 'zlib'

# The main class used to work with PNG file 
class PNG
  # PNG file's constant magic number
  SIGNATURE = "\x89PNG\r\n\x1A\n".force_encoding("BINARY")
  # Constant stop byte that will be added to the end of encoded sequence
  STOP_BYTE = "$"
  
  # Class attributes are extensive to make debugging easier
  attr_reader :filename, :user_bit_depth, :user_compression_level, :user_alpha, :datastream, :width, :height, :bit_depth, :colour_type, :compression_method, :filter_method, :interlace_method, :zlib_datastream, :zlib_compression_method, :zlib_compression_info, :zlib_check, :zlib_preset_dictionary, :zlib_compression_level, :zlib_checksum, :filtered_image_data, :image_data, :packed_image_data, :new_compressed_image_data, :verbose
  
  # Initialization of the PNG object
  def initialize(filename, user_bit_depth, user_compression_level, user_alpha, verbose)
    # Sets the level of verbose output requested by the interface module
    @verbose = verbose
    
    # Detailed information for verbose output
    info(verbose, 2, "Initializing PNG instance")
    
    # Sets the data provided by the interface
    @filename = filename
    @user_bit_depth = user_bit_depth
    @user_compression_level = user_compression_level
    @user_alpha = user_alpha

    # Opens the file
    file_open
    
    # Signature validation to verify that file is a PNG
    validate_signature

    # Reads the contents of the file
    read_file

    # Verifies that compressed image data is correct
    verify_zlib_datastream
  end

  # Decompresses the image data, verifies it and removes the filter, forming the raw data
  def unpack_image_data
    supported_file?
    decompress_image_data
    verify_filtered_image_data
    remove_filter
  end

  # Packs the image data, compresses it and writes to the file
  def save_file(output_file)
    pack_image_data
    compress_image_data
    write_file(output_file)
  end

  # Embeds data inside the image data
  def embed_data(data)
    # Verbose output
    info(verbose, 1, "Embedding data")

    # Adds the stop byte to the data string
    data = data + STOP_BYTE
    # Converts it to the string representing bits in the data
    data = data.unpack("B*").first

    # Loops through every pixel in the image data
    for y in 0..height-1
      for x in 0..width-1
        # Loops through every channel including alpha if requested by the user
        for channel in 0..(user_alpha ? 3 : 2)

          # Converts colour byte into a string representing its bits
          temp = image_data[y][x][channel].to_s(2).rjust(8, "0")

          # If the remaining data to write exceeds the affected bits depth
          if data.length > user_bit_depth
            # Moves bits from data string into the current colour byte
            temp[8-user_bit_depth..7] = data.slice!(0..user_bit_depth-1)
          # The very last few bits remains to be written
          else
            # Moves the last bits so no more data left after this
            temp[8-data.length..7] = data.slice!(0..data.length-1)
          end
          
          # Saves the edited colour byte into the image data
          image_data[y][x][channel] = temp.to_i(2)
          
          # Breaks from looping if no data left
          break if data.length.zero?
        end
        break if data.length.zero?
      end
      break if data.length.zero?
    end
  end

  # Extracts the data embedded into the image data till the stop byte
  def extract_data
    # Verbose output
    info(verbose, 1, "Extracting data")
    # Buffer is a string containing extracted bits
    buffer = ""
    # String containing extracted bytes
    extracted_data = ""
    # Flag indicating that stop byte was found
    stop_byte_found = false

    # Loops through every pixel in the image data
    for y in 0..height-1
      for x in 0..width-1
        # Loops through every channel including alpha if requested by the user
        for channel in 0..(user_alpha ? 3 : 2)

          # Saves bits extracted from the current colour byte to the buffer
          buffer << image_data[y][x][channel][0..user_bit_depth-1].to_s(2).rjust(user_bit_depth, "0")

          # In case if buffer contains enough bits to form a byte
          if buffer.length >= 8
            # Extracts the byte from the buffer
            extracted_byte = [buffer.slice!(0..7).to_i(2)].pack("C")

            # Checks if the extracted byte is the stop byte
            if extracted_byte == STOP_BYTE
              # Verbose output
              info(verbose, 2, "Stop byte found")
              # Sets the flag
              stop_byte_found = true 
            else
              # Otherwise just adds the extracted byte to the extracted data string
              extracted_data << extracted_byte
            end
          end

          # Breaks from looping if all of the data was extracted
          break if stop_byte_found
        end
        break if stop_byte_found
      end
      break if stop_byte_found
    end

    # If the stop byte was not found
    if !stop_byte_found
      # Tells the user that there is no hidden data
      puts "No data was found"
      exit 0
    end
    
    # Returns the extracted data
    extracted_data
  end
  
  # Verifies the size of the filtered image data
  def verify_filtered_image_data
    # Verbose output
    info(verbose, 2, "Verifying filtered image data")
    # Calculates the required size including filters
    expected_size = width * height * colour_channels + height
    # Error in case if the actual size of filtered image data is deifferent
    raise "Invalid filtered data size (#{filtered_image_data.length}), expected #{expected_size}" unless filtered_image_data.length == expected_size
  end

  # Verifies the filter according to the PNG specification
  def verify_filter(filter)
    raise "Invalid filter (#{filter})" unless filter < 5
  end
  
  # Function used by filter type 4 according to the PNG specification
  def paeth_predictor(a, b, c)
    p = a + b - c
    pa = (p - a).abs
    pb = (p - b).abs
    pc = (p - c).abs
    
    return a if pa <= pb and pa <= pc

    return b if pb <= pc

    c
  end

  # Saves the colour byte into image data after removing the filter
  def reconstruct(x, y, channel, filter)
    # Extracts the filtered colour byte from the filtered data
    filtered_byte = filtered_image_data[y*width*colour_channels + 1 + y + channel+x*colour_channels].unpack("C").first

    # Depending on the filter applies different reconstruction method
    case filter
    when 0
      # Filter type 0 does not modify the data
      reconstructed_byte = filtered_byte
    when 1
      # Pixel before the current one or 0 if current one is the most left pixel
      a = x.zero? ? 0 : image_data[y][x-1][channel]
      
      # Reconstructs the byte using previous pixel
      reconstructed_byte = filtered_byte + a
    when 2
      # Same as filter type 1 but uses pixel above instead
      b = y.zero? ? 0 : image_data[y-1][x][channel]
      reconstructed_byte = filtered_byte + b
    when 3
      a = x.zero? ? 0 : image_data[y][x-1][channel]
      b = y.zero? ? 0 : image_data[y-1][x][channel]
      # Adds the average of previous pixel and one above
      reconstructed_byte = filtered_byte + (a + b) / 2
    when 4
      a = x.zero? ? 0 : image_data[y][x-1][channel]
      b = y.zero? ? 0 : image_data[y-1][x][channel]
      # Pixel that is located diagonally to the top left
      c = (x.zero? or y.zero?) ? 0 : image_data[y-1][x-1][channel]
      # Relies on specific function
      reconstructed_byte = filtered_byte + paeth_predictor(a, b, c)
    end

    # Modulo 256 arithmetics prevents from overflowing
    reconstructed_byte = reconstructed_byte % 256

    # The reconstructed byte is saved into the image data
    image_data[y][x][channel] = reconstructed_byte
  end

  # Generates the raw image data
  def remove_filter
    # Verbose output
    info(verbose, 1, "Removing filters and reconstructing image data...")
    # Creates a 3D array that stores pixel data like this: [y][x][channel]
    @image_data = Array.new(height) { Array.new(width) { Array.new } }

    # Loops through every scanline
    for y in 0..height-1

      # First byte on each scanline indicates the filter applied to the whole scanline
      filter = filtered_image_data[y * (width * colour_channels + 1)].unpack("C").first

      # Verifies that filter type is valid
      verify_filter(filter)

      # Loops through every pixel in the current scanline
      for x in 0..width-1
        
        # Loops through every colour channel of the current pixel
        for channel in 0..colour_channels-1
          
          # Reconstructs the colour byte using filter from the scanline
          reconstruct(x, y, channel, filter)
        end
        
        # If the opened image does not contain alpha channel but user requested to use it
        if colour_channels == 3 and user_alpha
          # Adds additional channel indicating that current pixel iscompletely opaque
          image_data[y][x][3] = 255
        end
      end      
    end
    # If the opened image does not contain alpha channel but user requested to use it
    if colour_channels == 3 and user_alpha
      # Verbose output
      info(verbose, 2, "Alpha channel was added to the image data")
    end
  end

  # Packs that image data into the filtered image data string using filter type 0
  def pack_image_data
    # Verbose output
    info(verbose, 1, "Adding filters and packing image data...")

    # Empty string that will be used for storing filtered image data
    @packed_image_data = "".force_encoding("ASCII-8BIT")
    
    # Loops through each scanline in the image data
    for y in 0..height-1

      # Adds byte indicating filter type 0 for no modification of the raw data
      packed_image_data << 0#"\x00"#[0].pack("N").force_encoding("BINARY")

      # Loops through each pixel on the scanline
      for x in 0..width-1

        # Loops through every colour channel including alpha channel if requested
        for channel in 0..(user_alpha ? 3 : (colour_channels - 1))
         # Adds current byte from the image data to the filtered image data string
         packed_image_data << image_data[y][x][channel] 
        end
      end
    end
  end

  # Specifies Zlib FLEVEL compression value according to RFC-1950
  def get_compression_level
    case zlib_compression_level
    when 0
      "fastest algorithm"
    when 1
      "fast algorithm"
    when 2
      "default algorithm"
    when 3
      "maximum compression, slowest algorithm"
    end
  end
  
  # Verifies Zlib compression method according to RFC-1950
  def verify_zlib_compression_method
    # CM must be equal to 8
    raise "Invalid Zlib CM value (#{zlib_compression_method})" unless zlib_compression_method == 8
  end

  # Verifies Zlib compression info according to RFC-1950
  def verify_zlib_compression_info
    # CINFO must be less than 8
    raise "Invalid Zlib CINFO value (#{zlib_compression_info})" unless zlib_compression_info < 8
  end

  # Calculates LZ77 window size according to RFC-1950
  def get_lz77_window_size
    # Solution for CINFO = log2(LZ77_window_size) - 8
    2**(zlib_compression_info+8)
  end

  # Verifies Zlib check value according to RFC-1950
  def verify_zlib_check(cmf, flg)
    # CMF and FLG when viewed as 16-bit unsigned integer stored in MSB order must be a multiple of 31
    raise "Invalid Zlib FCHECK value (#{zlib_check})" unless ((cmf*256 + flg) % 31).zero?
  end

  # Verifies Zlib preset dictionary flag according to RFC-1950
  def verify_zlib_preset_dictionary
    # FDICT bit must be set to 0
    raise "Invalid Zlib FDICT flag (#{zlib_preset_dictionary})" unless !zlib_preset_dictionary
  end

  # Verifies Zlib compression level according to RFC-1950
  def verify_zlib_compression_level
    # FLEVEL must be less than 4
    raise "Invalid Zlib FLEVEL value (#{zlib_compression_level})" unless zlib_compression_level < 4
  end
  
  # Saves and verifies Zlib datastream parameters
  def verify_zlib_datastream
    # Verbose output
    info(verbose, 2, "Verifying Zlib datastream")

    # Offset for 32-bit checksum should be 16-bit header + compresssed data
    adler32_offset = zlib_datastream.length - 6

    # Unpacks CMF, FLG and ADLER32 according to RFC-1950
    cmf, flg, @zlib_checksum = zlib_datastream.unpack("CCx#{adler32_offset}N")

    # Extracts CM from CMF
    @zlib_compression_method = cmf[0..3]
    verify_zlib_compression_method

    # Extracts CINFO from CMF
    @zlib_compression_info = cmf[4..7]
    verify_zlib_compression_info

    # Extracts FCHECK value from FLG
    @zlib_check = flg[0..4]
    verify_zlib_check(cmf, flg)

    # Extracts FDICT flag from FLG
    @zlib_preset_dictionary = !flg[5].zero?
    verify_zlib_preset_dictionary

    # Extracts FLEVEL value from FLG
    @zlib_compression_level = flg[6..7]
    verify_zlib_compression_level
  end

  # Decompresses filtered image data
  def decompress_image_data
    # Verbose output
    info(verbose, 1, "Decompressing image data")

    # Inflates Zlib datastream, i.e. decompresses filtered image data
    @filtered_image_data = Zlib.inflate(zlib_datastream)

    # Verifies its integrity
    verify_decompressed_image_data
  end
  
  # Verifies the integrity of the inflated filtered image data
  def verify_decompressed_image_data
    # Verbose output
    info(verbose, 2, "Verifying decompressed image data")

    # Raises error if calculated checksum is not equal to the provided
    raise "Adler-32 mismatch" unless @zlib_checksum == Zlib.adler32(@filtered_image_data)
  end

  # Compresses packed (filtered) image data
  def compress_image_data
    # Verbose output
    info(verbose, 1, "Compressing image data")

    # Deflates packed image data using requested compression level (0 - min, 9 - max)
    @new_compressed_image_data = Zlib.deflate(packed_image_data, user_compression_level)
  end

  # Returns the number of colour channels in the opened image
  # Might be 3 even if user requests to use alpha channel
  def colour_channels
    # For truecolour it is either 3 (RGB) or 4 (RGBA)
    colour_type == 2 ? 3 : 4
  end
  
  # Saves (modified) PNG data into the file
  def write_file(output_file)
    # Verbose output
    info(verbose, 1, "Saving new image to the output file")

    # Creates new file
    outfile = File.new(output_file, "wb")
    # Adds PNG signature
    outfile.write(SIGNATURE)
    # Writes IHDR, IDAT and IEND chunks
    outfile.write(write_chunk("IHDR")+write_chunk("IDAT", new_compressed_image_data)+write_chunk("IEND"))
    outfile.close
  end

  # Calculates CRC checksum for the chunk of type "type" with content "content"
  def calculate_crc(type, content)
    # Relies on Zlib's method
    Zlib.crc32(type + content)
  end

  # Returns string that is a valid chunk of type "type" optionally containing "content"
  def write_chunk(type, content = "")
    # Verbose output
    info(verbose, 2, "Writing #{type} chunk")
    # For chunk IHDR it will create the content from the stored information
    if type == "IHDR"
      # Packs the IHDR data modifying colour type if alpha channel was forced onto RGB
      content = [width, height, bit_depth, (user_alpha ? 6 : (colour_type == 6 ? 6 : 2)), compression_method, filter_method, interlace_method].pack("NNC5")
    end
    # Calculates the length of the content
    length = [content.length].pack("N")
    # Calculates the CRC checksum for the chunk
    crc = [calculate_crc(type, content)].pack("N")

    # Returns the string in the proper chunk format
    length + type + content + crc
  end

  # Prints information about the provided PNG file telling how much data it can store
  def print_info
    puts "File: \"#{filename}\""
    puts "-- PNG Parameters --"
    puts "Dimensions: #{width}x#{height}"
    puts "Bit Depth: #{bit_depth}"
    puts "Colour Type: #{colour_type} (#{get_image_type})"
    puts "Compression Method: #{compression_method}"
    puts "Filter Method: #{filter_method}"
    puts "Interlace Method: #{interlace_method}"
    puts "-- Zlib Parameters --"
    puts "Compression Method: #{zlib_compression_method} (Deflate)"
    puts "Compression Info: #{zlib_compression_info} (LZ77 window size is #{get_lz77_window_size} bytes)"
    puts "Check Value: #{zlib_check}"
    puts "Preset Dictionary: #{zlib_preset_dictionary}"
    puts "Compression Level: #{zlib_compression_level} (#{get_compression_level})"
    puts "-- User Parameters --"
    puts "Colour Channels: RGB" + (user_alpha ? "A" : "")
    puts "Bit Depth: #{user_bit_depth}"
    puts "Possible to hide #{available_size} bytes"
  end

  # Calculates size in bytes of how much data could be embedded into the file
  def available_size
    # Size is calculated considering request for bit depth and alpha channel
    size = width*height*(user_alpha ? 4 : 3)*user_bit_depth/8
    # Minus 1 byte for the stop byte
    size = size - 1 if size > 0
    size
  end
  
  # Returns the image type name according to the PNG specification
  def get_image_type
    case colour_type
    when 0
      "Greyscale"
    when 2
      "Truecolour"
    when 3
      "Indexed-colour"
    when 4
      "Greyscale with alpha"
    when 6
      "Truecolour with alpha"
    end
  end

  # Checks if PNG file is supported by the MeiSai
  def supported_file?
    # Verbose output
    info(verbose, 2, "Checking if PNG file is supported")
    # Raises the error if PNG file is not 8 bit depth or is not using truetype colour mode
    raise "Expecting 8-bit colour depth" unless bit_depth == 8
    raise "Expecting truecolour with or without alpha" unless colour_type == 2 || colour_type == 6
    # Additionally checks that filter and compression methods are correct
    raise "Expecting filter method 0" unless filter_method == 0
    raise "Expecting deflate compression (method 0)" unless compression_method == 0
  end

  # Opens the file creating IO datastream
  def file_open
    # Verbose output
    info(verbose, 1, "Openning PNG file")
    @datastream = File.open(filename, "rb")
  end  
  
  # Validates the file's PNG signature
  def validate_signature
    # Verbose output
    info(verbose, 2, "Verifying PNG signature")
    # Reads the signature from the file
    signature = datastream.read(SIGNATURE.length)
    
    raise "Invalid PNG signature" unless signature == SIGNATURE
  end

  # Parses the remaining contents of the file in the datastream (without the signature)
  def read_file
    # Verbose output
    info(verbose, 1, "Reading contents of the file")

    # Creates a string to store contents of the IDAT chunk(s)
    # It will contain Zlib datastream that is a compressed (deflated) filtered image data
    @zlib_datastream = ""

    # Reads file chunk by chunk till the end of the datastream
    until datastream.eof?
      # Invokes the method to read current chunk
      type, content = read_chunk
      
      # Depending on the chunk type performs specific actions
      case type
      # Header chunk containing information about the image data
      when "IHDR"
        # IHDR chunk's content length is always 13 Bytes, must be parsed in certain manner
        properties = content.unpack("NNC5")
        # These properties are saved for later
        @width, @height, @bit_depth, @colour_type, @compression_method, @filter_method, @interlace_method = properties
      # Chunk containing compressed image data
      when "IDAT"
        # Verbose output
        info(verbose, 1, "Extracting compressed image data")
        # Might be multiple IDAT chunks, so appends the contents to the Zlib's datastream
        @zlib_datastream << content
      end
    end
  end

  # Reads the chunk
  def read_chunk
    # First 4 bytes indicate the length of the chunk's content, next 4 bytes - the name
    length, type = datastream.read(8).unpack("NA4")
    # Verbose output
    info(verbose, 2, "Reading #{type} chunk")
    # Next is the chunk's content of the specified length
    content = datastream.read(length)
    
    # Content is followed by 4 byte checksum used to validate the integrity of the chunk
    validate_chunk(type, content)

    # Returns the type and content
    [type, content]
  end

  # Validates the chunk from datastream
  def validate_chunk(type, content)
    # Verbose output
    info(verbose, 2, "Validating #{type} chunk")
    # Reads the 32-bit CRC checksum
    crc_read = datastream.read(4).unpack("N").first
    
    # Calculates the checksum for the chunk
    crc_calculated = calculate_crc(type, content)
    
    # Raises the error if they mismatch
    raise "CRC mismatch for chunk #{type}" unless crc_read == crc_calculated
  end
end