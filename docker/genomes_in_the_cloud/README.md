
# Login to ECR
```bash
$(aws ecr get-login --no-include-email --region us-east-1)
```

# Build the image and push it to the repository
```bash
export REPO=123456789.dkr.ecr.us-east-1.amazonaws.com/genomes-in-the-cloud
export VERSION=1.0.0

docker build -t ${REPO}:${VERSION} --rm .
docker push ${REPO}:${VERSION}
```
