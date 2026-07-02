# Differential-harness image: compiles the C LaSRC leaves against the conda
# espa-product-formatter / hdf headers and runs the C-vs-Rust proptests.
#
# Mirrors hls-science-container's pixi-based build, minus the Fmask/MCR bulk.
# Build context is the REPO ROOT (needs lasrc/c_version/src + lasrc_core +
# lasrc_cref), Dockerfile lives in lasrc_cref/.
#
#   docker build -f lasrc_cref/Dockerfile --platform linux/amd64 \
#     --secret id=aws_access_key_id,env=AWS_ACCESS_KEY_ID \
#     --secret id=aws_secret_access_key,env=AWS_SECRET_ACCESS_KEY \
#     --secret id=aws_session_token,env=AWS_SESSION_TOKEN \
#     -t lasrc-cref .
#   docker run --rm lasrc-cref            # runs `pixi run test-lut`
#
# syntax=docker/dockerfile:1.7-labs
ARG PLATFORM=linux/amd64
FROM --platform=${PLATFORM} ghcr.io/prefix-dev/pixi AS build

WORKDIR /app

# Only the pieces the harness needs. Paths match the crate's relative refs:
# build.rs -> ../lasrc/c_version/src, Cargo.toml -> ../lasrc_core.
COPY lasrc/c_version/src ./lasrc/c_version/src
COPY lasrc_core ./lasrc_core
COPY lasrc_cref ./lasrc_cref

WORKDIR /app/lasrc_cref

# espa-product-formatter + espa-surface-reflectance live on the private
# s3://hls-conda-channels channel; pixi needs AWS creds to fetch them.
RUN --mount=type=cache,target=/root/.cache/rattler/cache \
    --mount=type=secret,id=aws_access_key_id,required=false \
    --mount=type=secret,id=aws_secret_access_key,required=false \
    --mount=type=secret,id=aws_session_token,required=false \
    AWS_DEFAULT_REGION=us-west-2 AWS_REGION=us-west-2 \
    AWS_ACCESS_KEY_ID=$(cat /run/secrets/aws_access_key_id 2>/dev/null) \
    AWS_SECRET_ACCESS_KEY=$(cat /run/secrets/aws_secret_access_key 2>/dev/null) \
    AWS_SESSION_TOKEN=$(cat /run/secrets/aws_session_token 2>/dev/null) \
    pixi install

# Stage 1+ differential tests. The task sets PREFIX=$CONDA_PREFIX so cc finds
# the conda headers/libs, and runs inside the activated env so it inherits the
# same $CC/$CFLAGS used for the shipped binary.
CMD ["pixi", "run", "test-lut"]
