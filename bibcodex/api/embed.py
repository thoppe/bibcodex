from typing import List, Dict
import numpy as np
import torch
from .helpers import CachedDownloader


class HF_transformer(CachedDownloader):
    """
    Uses huggingface transformers to run an embeding of the text.
    """

    name = "embed"
    chunksize = 10
    datatype = np.ndarray

    def __init__(self):

        self.model = None
        self.tokenizer = None
        self.model_name = "allenai/specter"

        # Set device on GPU if available, else CPU
        self.device = torch.device(
            "cuda" if torch.cuda.is_available() else "cpu"
        )

        super().__init__()

    def _load_model(self):
        # Lazy import of the transformer model since it's expensive

        if self.model is None:
            from transformers import AutoTokenizer, AutoModel

            name = self.model_name

            self.tokenizer = AutoTokenizer.from_pretrained(name)
            self.model = AutoModel.from_pretrained(name)

            # Put the model in evaluation mode and move to device
            self.model.eval().to(self.device)

    @CachedDownloader.cached
    def get_from_text(
        self,
        text: List[str],
    ) -> Dict[str, datatype]:

        self._load_model()

        tokens = self.tokenizer(
            text,
            padding=True,
            truncation=True,
            return_tensors="pt",
            max_length=512,
        ).to(self.device)

        result = self.model(**tokens)

        # Take the first token in the batch as the embedding
        embeddings = result["last_hidden_state"][:, 0, :]

        # Convert the embeddings from torch to a numpy array
        embeddings = embeddings.detach().cpu().numpy()

        # Zip the data back up for the caching
        data = dict(zip(text, embeddings))

        return data


downloader = HF_transformer()
