from fastapi import FastAPI
from app.routers import sequences, annotations
from app.core.config import settings

app = FastAPI(
    title="BioSequence API",
    description="A bioinformatics API for DNA/RNA sequence analysis and annotation identification",
    version="1.0.0",
)

# Include routers
app.include_router(sequences.router, prefix="/api/v1/sequences", tags=["sequences"])
app.include_router(
    annotations.router, prefix="/api/v1/annotations", tags=["annotations"]
)


@app.get("/")
async def root():
    return {"message": "Welcome to BioSequence API", "version": "1.0.0"}


@app.get("/health")
async def health_check():
    return {"status": "healthy"}


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)
