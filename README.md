# Square Detection Assignment

This project processes images to detect and work with squares in grayscale images. The assignment is divided into three main sections:

## Section 1: Find Squares
- **Input**: An image containing multiple non-overlapping squares with grayscale values between 1-255. The background is 0.
- **Goal**: Traverse the image and return a list of detected squares with the following properties:
  - Top-left corner coordinates.
  - Edge width.
  - Grayscale value.
  
## Section 2: Recursive Square Detection
- **Goal**: Repeat the detection from Section 1, but implement a **recursive** method for finding squares. This approach uses DFS, which may cause performance issues for large images.

## Section 3: Diagonal Square Overlap
- **Input**: An image with a diagonal square and a larger image with horizontal squares.
- **Goal**: Detect the diagonal square and find the horizontal square with the best overlapping fit, based on maximum area overlap.

## Square Class
The `Square` class represents a square with the following properties:
- Row and column of the top-left corner.
- Edge size.
- Grayscale value.

### Example Usage:
```cpp
vector<Square> squares = findAllSquares(image);
Square* bestFit = findBestOverlap(diagonalSquare, squares);
