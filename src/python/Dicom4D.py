import struct
import numpy as np
import nibabel as nib

class Dicom4D:
    """
    """
    def __init__ (self, fndcm):
        self.Filename = fndcm

        # Read file into buffer
        self.Buffer = Dicom4DBuffer()
        self.__load()
        
    # Load buffer
    def __load (self):
        buf = self.Buffer

        # Get Header Info
        f = open(self.Filename, 'rb')
        
        # Skip first 128 Bytes
        f.seek(128)
        
        # DICOM header
        DICM = f.read(4).decode('utf-8')
        print("Header: ", DICM)

        # The data tag (7fe0, 0010)
        dataTag = (0x7fe0, 0x0010)
        
        LoopCnt = 0

        # Initialize current tag
        tag = (0x0000, 0x0000)

        while (tag != dataTag and LoopCnt < 200):
            LoopCnt += 1

            tag = (int.from_bytes(f.read(2), byteorder = 'little'), int.from_bytes(f.read(2), byteorder = 'little'))

            code = f.read(2).decode('utf-8')
            n = int.from_bytes(f.read(2), byteorder = 'little')
            # print('(', '{:04x}'.format(tag[0]), '{:04x}'.format(tag[1]), ')', 'code = ', code, '\t n = ', n)
            
            if (tag[0] == 0x0018):
                if (tag[1] == 0x602c):
                    buf.deltaX = struct.unpack('<d', f.read(n))[0] * 10
                    # print('deltaX = ', buf.deltaX)
                elif (tag[1] == 0x602e):
                    buf.deltaY = struct.unpack('<d', f.read(n))[0] * 10
                    # print('deltaY = ', buf.deltaY)
                elif (tag[1] == 0x1063):
                    # print('frameTime size = ', n)
                    buf.frameTime = float(f.read(n).decode('utf-8'))
                    # print('frameTime = ', buf.frameTime)
                else:
                    f.read(n)
            elif (tag[0] == 0x0028):
                if (tag[1] == 0x0008):
                    buf.numFrames = int(f.read(n).decode('utf-8'))
                    # print('numFrames = ', buf.numFrames)
                elif (tag[1] == 0x0010):
                    buf.height = int.from_bytes(f.read(n), byteorder = 'little')
                    # print('height = ', buf.height)
                elif (tag[1] == 0x0011):
                    buf.width = int.from_bytes(f.read(n), byteorder = 'little')
                    # print('width = ', buf.width)
                else:
                    f.read(n)
            elif (tag[0] == 0x3001):
                if (tag[1] == 0x1001):
                    buf.depth = int.from_bytes(f.read(n), byteorder = 'little')
                    # print('depth = ', buf.depth)
                elif (tag[1] == 0x1003):
                    buf.deltaZ = struct.unpack('<d', f.read(n))[0] * 10
                    # print('deltaZ = ', buf.deltaZ)
                else:
                    f.read(n)
            else:
                f.read(n)

            if (code == 'OB'):
                f.read(6)

        # Load Data
        n = int.from_bytes(f.read(4), byteorder = 'little')
        expectedSize = buf.width * buf.height * buf.depth * buf.numFrames
        
        # Data length should match total # of voxels
        assert(n == expectedSize)

        raw = f.read(n)
        assert(len(raw) == expectedSize)
        
        # Import buffer into a numpy array
        # We need to reverse the dimension order (t*z*y*x) to load the 1D buffer,
        # and transpose the array (x*y*z*t) to adjust to the order that numpy 
        # organizes array dimensions
        bufferArr = np.frombuffer(raw, np.uint8) \
            .reshape((buf.numFrames, buf.depth, buf.height, buf.width))
        bufferArr = np.transpose(bufferArr)

        buf.data = bufferArr
        buf.printInfo()

        f.close()

    # Export all frames from buffer
    def Export4D (self, outfn):
        print("Exporting 4D Image to: ", outfn)
        buf = self.Buffer
        affine = buf.GetAffine()
        img = nib.Nifti1Image(buf.data, affine)
        # print(img)

        # Save image to the file
        nib.save(img, outfn)


    # Export specific frame
    def ExportFrame (self, frameNum, outfn):
        print("Exporting Frame ", frameNum, " to: ", outfn)
        buf = self.Buffer
        affine = buf.GetAffine()
        dataArr = np.transpose(np.transpose(buf.data)[frameNum])
        img = nib.Nifti1Image(dataArr, affine)

        # Save image to the file
        nib.save(img, outfn)


    def printInfo (self):
        print("Class: Dicom4DReader")
        print("Filename: ", self.Filename)


class Dicom4DBuffer:
    """
    """
    def __init__ (self):
        # Pixel Data
        self.data       = None
        # Dimension
        self.width      = None
        self.height     = None
        self.depth      = None
        self.numFrames  = None
        # Spacing
        self.deltaX     = None
        self.deltaY     = None
        self.deltaZ     = None
        self.frameTime  = None

    def GetAffine (self):
        return np.diag([self.deltaX, self.deltaY, self.deltaZ, 1])

    def printInfo (self):
        print("Class: Dicom4DData")
        print("Dimension: ", self.width, "x", self.height, "x", self.depth, "x", self.numFrames)
        print("Spacing: [", self.deltaX, ", ", self.deltaY, ", ", self.deltaZ, ", ", self.frameTime, "]")

    