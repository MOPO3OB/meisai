# This is the Base64 module (part of MeiSai Framework)

# Loads system's base64 library
require "base64"

# Encoding function
def encode(data, verbose)
  info(verbose, 1, "Encoding data")
  Base64.strict_encode64(data)
end

# Decoding function
def decode(data, verbose)
  begin
    info(verbose, 1, "Decoding data")
    # Returns the decoded data
    Base64.strict_decode64(data)
  # Catches error indicating failed attempt to decode the data
  rescue ArgumentError => e
    # Notifies user even on quite mode
    info(verbose, 0, "Could not decode data")
    exit 0
  end
end
